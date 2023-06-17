import os
import subprocess
import time

import click
import pathops
from cffsubr import subroutinize
from fontTools.fontBuilder import FontBuilder
from fontTools.misc.cliTools import makeOutputFileName
from fontTools.pens.qu2cuPen import Qu2CuPen
from fontTools.pens.recordingPen import DecomposingRecordingPen
from fontTools.pens.t2CharStringPen import T2CharStringPen
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools.subset import Subsetter
from fontTools.ttLib.scaleUpem import scale_upem
from fontTools.ttLib.ttFont import TTFont, TTLibError


class TrueTypeToCFFOptions(object):
    def __init__(self):
        self.tolerance: float = 1.0
        self.subroutinize = True
        self.remove_glyphs = True
        self.new_upem = None
        self.recalc_timestamp = False
        self.output_dir = None
        self.overwrite = True


class TrueTypeToCFFRunner(object):
    def __init__(self, fonts: list[TTFont]):
        super().__init__()
        self.fonts = fonts
        self.options = TrueTypeToCFFOptions()

    def run(self) -> None:
        converted_files_count = 0
        start_time = time.time()

        for count, source_font in enumerate(self.fonts, start=1):
            t = time.time()

            try:
                print()
                generic_info_message(f"Converting file {count} of {len(self.fonts)}: {source_font.reader.file.name}")

                ext = ".otf" if source_font.flavor is None else "." + str(source_font.flavor)
                # Suffix is necessary because without it source woff and woff2 files would be overwritten
                suffix = "" if source_font.flavor is None else ".otf"
                output_file = makeOutputFileName(
                    source_font.reader.file.name,
                    suffix=suffix,
                    extension=ext,
                    outputDir=self.options.output_dir,
                    overWrite=self.options.overwrite,
                )

                tolerance = self.options.tolerance / 1000 * source_font["head"].unitsPerEm
                source_font.recalcTimestamp = self.options.recalc_timestamp

                # Prepare the font for conversion
                decomponentize(font=source_font)
                if self.options.new_upem is not None:
                    scale_upem(source_font, new_upem=self.options.new_upem)
                if self.options.remove_glyphs:
                    remove_unneeded_glyphs(source_font)

                charstrings, glyphs_with_errors = get_qu2cu_charstrings_without_pathops(font=source_font, tolerance=tolerance)

                # Fix charstrings in case of errors
                if len(glyphs_with_errors) > 0:
                    generic_warning_message(f"The following glyphs must be fixed: {', '.join(glyphs_with_errors)}")

                    # Save the source font to a temp source file in case scale_upm and/or remove_unneeded_glyphs are
                    # True
                    temp_source_file = makeOutputFileName(output_file, suffix="_tmp", extension=".ttf", overWrite=True)
                    source_font.save(temp_source_file)

                    # Dump the temp source file as .cff
                    temp_cff_file = makeOutputFileName(output_file, suffix="_tmp", extension=".cff", overWrite=True)
                    run_shell_command(
                        args=["tx", "-cff", "-n", "+b", "-S", temp_source_file, temp_cff_file], suppress_output=True
                    )

                    # Build a temporary otf font
                    temp_otf_font = build_cff_font(source_font, charstrings)
                    temp_otf_file = makeOutputFileName(output_file, suffix="_tmp", extension=".otf", overWrite=True)
                    temp_otf_font.save(temp_otf_file)

                    # Merge the CFF dumped before into the temp otf font
                    run_shell_command(args=["sfntedit", "-a", f"CFF={temp_cff_file}", temp_otf_file])

                    # Get the charstrings from the temp otf font
                    temp_otf_font_1 = TTFont(temp_otf_file)
                    temp_charstrings, glyphs_with_errors_2 = get_t2_charstrings_with_pathops(temp_otf_font_1)

                    # Replace the charstrings with errors with the ones from the temporary otf font
                    for g in glyphs_with_errors:
                        try:
                            charstrings[g] = temp_charstrings[g]
                        except KeyError:
                            glyphs_with_errors_2.append(g)

                    # Display an error message if not all glyphs have been fixed
                    if len(glyphs_with_errors_2) > 0:
                        generic_error_message(
                            f"The following glyphs couldn't be fixed: {', '.join(glyphs_with_errors_2)}"
                        )
                    else:
                        generic_info_message("All glyphs have been fixed")

                    # Remove the temporary files
                    os.remove(temp_source_file)
                    os.remove(temp_otf_file)
                    os.remove(temp_cff_file)

                cff_font = build_cff_font(font=source_font, charstrings=charstrings)

                if self.options.subroutinize:
                    # cffsubr doesn't work with woff/woff2 fonts, we need to temporary set flavor to None
                    flavor = cff_font.flavor
                    if flavor is not None:
                        cff_font.flavor = None
                    subroutinize(cff_font)
                    cff_font.flavor = flavor

                cff_font.save(output_file)

                generic_info_message(f"Elapsed time: {round(time.time() - t, 3)} seconds")
                file_saved_message(output_file)
                converted_files_count += 1

            except Exception as e:
                generic_error_message(e)

        print()
        generic_info_message(f"Total files     : {len(self.fonts)}")
        generic_info_message(f"Converted files : {converted_files_count}")
        generic_info_message(f"Elapsed time    : {round(time.time() - start_time, 3)} seconds")


def build_cff_font(font: TTFont, charstrings: dict) -> TTFont:
    cff_font_info = get_cff_font_info(font)
    post_values = get_post_values(font)

    fb = FontBuilder(font=font)
    fb.isTTF = False
    for table in ["glyf", "cvt ", "loca", "fpgm", "prep", "gasp", "LTSH", "hdmx"]:
        if table in fb.font:
            del fb.font[table]
    fb.setupCFF(
        psName=font["name"].getDebugName(6),
        charStringsDict=charstrings,
        fontInfo=cff_font_info,
        privateDict={},
    )

    fb.setupDummyDSIG()
    fb.setupMaxp()
    fb.setupPost(**post_values)

    return fb.font


def get_qu2cu_charstrings_without_pathops(font: TTFont, tolerance: float = 1.0):
    charstrings = {}
    glyphs_with_errors = []
    glyph_set = font.getGlyphSet()

    for k, v in glyph_set.items():
        try:
            t2_pen = T2CharStringPen(v.width, glyphSet=glyph_set)
            qu2cu_pen = Qu2CuPen(t2_pen, max_err=tolerance, all_cubic=True, reverse_direction=True)
            glyph_set[k].draw(qu2cu_pen)
            charstring = t2_pen.getCharString()
            charstrings[k] = charstring
        except Exception as e:
            glyphs_with_errors.append(k)
            generic_warning_message(f"{k}: {e}")

    return charstrings, glyphs_with_errors


def get_qu2cu_charstrings(font: TTFont, tolerance: float = 1.0) -> (dict, list):
    charstrings = {}
    glyph_set = font.getGlyphSet()
    glyphs_with_errors = []

    for k, v in glyph_set.items():
        pathops_path = pathops.Path()
        pathops_pen = pathops_path.getPen(glyphSet=glyph_set)
        t2_pen = T2CharStringPen(v.width, glyphSet=glyph_set)
        qu2cu_pen = Qu2CuPen(t2_pen, max_err=tolerance, all_cubic=True, reverse_direction=False)

        try:
            glyph_set[k].draw(pathops_pen)
            pathops_path.simplify()
            pathops_path.draw(qu2cu_pen)

        except TypeError as e:
            generic_warning_message(f"{e} [{k}]")
            glyphs_with_errors.append(k)

        except NotImplementedError as e:
            generic_warning_message(f"{e}: [{k}]")
            glyphs_with_errors.append(k)
            qu2cu_pen.all_cubic = False
            glyph_set[k].draw(qu2cu_pen)

        charstring = t2_pen.getCharString()
        charstrings[k] = charstring

    return charstrings, glyphs_with_errors


def get_t2_charstrings(font: TTFont) -> dict:
    """
    Get CFF charstrings without pathops
    """
    charstrings = {}
    glyph_set = font.getGlyphSet()
    for k, v in glyph_set.items():
        t2_pen = T2CharStringPen(v.width, glyphSet=glyph_set)
        glyph_set[k].draw(t2_pen)
        charstrings[k] = t2_pen.getCharString()
    return charstrings


def get_t2_charstrings_with_pathops(font: TTFont) -> (dict, list):
    """
    Get CFF charstrings using T2CharStringPen and pathops

    :return: CFF charstrings.
    """
    charstrings = {}
    glyph_set = font.getGlyphSet()
    glyphs_with_errors = []

    for k, v in glyph_set.items():
        pathops_path = pathops.Path()
        pathops_pen = pathops_path.getPen(glyphSet=glyph_set)
        t2_pen = T2CharStringPen(v.width, glyphSet=glyph_set)
        try:
            glyph_set[k].draw(pathops_pen)
            pathops_path.simplify()
            pathops_path.draw(t2_pen)
        except TypeError:
            glyphs_with_errors.append(k)
            glyph_set[k].draw(t2_pen)

        charstring = t2_pen.getCharString()
        charstrings[k] = charstring

    return charstrings, glyphs_with_errors


def get_cff_font_info(font: TTFont) -> dict:
    """
    Setup CFF topDict

    :return: A dictionary of the font info.
    """

    font_revision = str(round(font["head"].fontRevision, 3)).split(".")
    major_version = str(font_revision[0])
    minor_version = str(font_revision[1]).ljust(3, "0")

    cff_font_info = dict(
        version=".".join([major_version, str(int(minor_version))]),
        FullName=font["name"].getBestFullName(),
        FamilyName=font["name"].getBestFamilyName(),
        ItalicAngle=font["post"].italicAngle,
        UnderlinePosition=font["post"].underlinePosition,
        UnderlineThickness=font["post"].underlineThickness,
        isFixedPitch=False if font["post"].isFixedPitch == 0 else True,
    )

    return cff_font_info


def get_post_values(font: TTFont) -> dict:
    post_info = dict(
        italicAngle=round(font["post"].italicAngle),
        underlinePosition=font["post"].underlinePosition,
        underlineThickness=font["post"].underlineThickness,
        isFixedPitch=font["post"].isFixedPitch,
        minMemType42=font["post"].minMemType42,
        maxMemType42=font["post"].maxMemType42,
        minMemType1=font["post"].minMemType1,
        maxMemType1=font["post"].maxMemType1,
    )
    return post_info


def run_shell_command(args, suppress_output=False):
    """
    Runs a shell command.
    Returns True if the command was successful, and False otherwise.
    """
    sup = subprocess.DEVNULL if suppress_output else None

    try:
        subprocess.check_call(args, stderr=sup, stdout=sup)
        return True
    except (subprocess.CalledProcessError, OSError) as err:
        print(err)
        return False


def decomponentize(font: TTFont) -> None:
    if "glyf" not in font:
        return

    glyph_set = font.getGlyphSet()
    glyf_table = font["glyf"]
    dr_pen = DecomposingRecordingPen(glyph_set)
    tt_pen = TTGlyphPen(None)

    for glyph_name in font.glyphOrder:
        glyph = glyf_table[glyph_name]
        if not glyph.isComposite():
            continue

        dr_pen.value = []
        tt_pen.init()

        glyph.draw(dr_pen, glyf_table)
        dr_pen.replay(tt_pen)

        glyf_table[glyph_name] = tt_pen.glyph()


def remove_unneeded_glyphs(font: TTFont) -> None:
    glyph_ids_to_remove = []
    for g in [".null", "NUL", "NULL", "uni0000", "CR", "nonmarkingreturn", "uni000D"]:
        try:
            glyph_ids_to_remove.append(font.getGlyphID(g))
        except KeyError:
            pass

    glyph_ids = [i for i in font.getReverseGlyphMap().values() if i not in glyph_ids_to_remove]
    if len(glyph_ids_to_remove) > 0:
        subsetter = Subsetter()
        subsetter.options.drop_tables = []
        subsetter.options.passthrough_tables = True
        subsetter.options.name_IDs = "*"
        subsetter.options.name_legacy = True
        subsetter.options.name_languages = "*"
        subsetter.options.layout_features = "*"
        subsetter.options.hinting = False
        subsetter.options.notdef_glyph = True
        subsetter.options.notdef_outline = True
        subsetter.glyph_ids_requested = glyph_ids
        subsetter.subset(font)


def get_fonts_in_path(input_path: str) -> list[TTFont]:
    """
    Returns a list of TTFont objects in input_path

    :param input_path: path to file or folder
    :return: a list of TTFont objects
    """
    files = []
    if os.path.isfile(input_path):
        files.append(input_path)
    if os.path.isdir(input_path):
        files = [os.path.join(input_path, os.path.basename(f)) for f in os.listdir(input_path)]

    fonts = []
    for file in files:
        try:
            font = TTFont(file)
            # Filter only static TrueType fonts
            is_true_type = "glyf" in font
            is_static = "fvar" not in font
            if is_true_type and is_static:
                fonts.append(font)
        except (TTLibError, PermissionError, Exception):
            pass

    return fonts


def get_output_dir(input_path: str, output_dir: str = None) -> str:
    """
    If the output directory is not specified, then the output directory is the directory of the input file if the input
    is a file, or the input directory if the input is a directory

    :param input_path: The path to the input file or directory
    :type input_path: str
    :param output_dir: The output directory, if specified
    :type output_dir: str
    :return: The output directory.
    """
    if output_dir is not None:
        return os.path.abspath(output_dir)
    else:
        if os.path.isfile(input_path):
            return os.path.dirname(input_path)
        else:
            return input_path


def file_saved_message(file: str):
    click.secho(f"[{click.style('DONE', fg='green')}] {file} {click.style('saved', fg='green')}")


def generic_error_message(error_message):
    click.secho(f"[{click.style('FAIL', fg='red')}] {error_message}")


def generic_info_message(info_message, nl=True):
    click.secho(f"[{click.style('INFO', fg='cyan')}] {info_message}", nl=nl)


def generic_warning_message(info_message, nl=True):
    click.secho(f"[{click.style('WARN', fg='yellow')}] {info_message}", nl=nl)


@click.command()
@click.argument(
    "input_path",
    type=click.Path(exists=True, resolve_path=True, dir_okay=True, file_okay=True),
)
@click.option(
    "-t",
    "--tolerance",
    type=click.FloatRange(0, 3),
    default=1,
    help="""
    Conversion tolerance (0-3, default 1). Low tolerance adds more points but keeps shapes. High tolerance adds few
    points but may change shape.
    """,
)
@click.option(
    "--scale-upm",
    "new_upem",
    type=click.IntRange(1000, 8192, max_open=True),
    default=None,
    help="""
    Scale units-per-em of source font to passed value before converting to OTF
    """,
)
@click.option(
    "--remove-glyphs",
    is_flag=True,
    help="""
    Remove NULL and CR glyphs from output fonts
    """,
)
@click.option(
    "--no-subr",
    "apply_subroutines",
    is_flag=True,
    default=True,
    help="""
    Do not subroutinize converted fonts with cffsubr
    """,
)
@click.option(
    "-out",
    "--output-dir",
    "outputDir",
    type=click.Path(file_okay=False, resolve_path=True),
    default=None,
    help="""
    Specify the directory where output files are to be saved. If output_dir doesn't exist, will be created. If not
    specified, files are saved to the source folder.
    """,
)
@click.option(
    "--recalc-timestamp",
    "recalcTimestamp",
    is_flag=True,
    default=False,
    help="""
    Keep the original font 'modified' timestamp (head.modified) or set it to current time. By default, original
    timestamp is kept.
    """,
)
@click.option(
    "--no-overwrite",
    "overWrite",
    is_flag=True,
    default=True,
    help="""
    Overwrite existing output files or save them to a new file (numbers are appended at the end of file name). By
    default, files are overwritten.
    """,
)
def cli(
        input_path,
        tolerance=1,
        new_upem=None,
        remove_glyphs=False,
        apply_subroutines=True,
        outputDir=None,
        recalcTimestamp=False,
        overWrite=True,
):
    """
    Converts TTF fonts (or TrueType flavored woff/woff2 web fonts) to OTF fonts (or CFF flavored woff/woff2 web fonts).
    """

    fonts = get_fonts_in_path(input_path)
    output_dir = get_output_dir(input_path=input_path, output_dir=outputDir)

    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception as e:
        generic_error_message(e)
        return

    converter = TrueTypeToCFFRunner(fonts)
    converter.options.tolerance = tolerance
    converter.options.new_upem = new_upem
    converter.options.remove_glyphs = remove_glyphs
    converter.options.subroutinize = apply_subroutines
    converter.options.recalc_timestamp = recalcTimestamp
    converter.options.output_dir = output_dir
    converter.options.overwrite = overWrite
    converter.run()


if __name__ == "__main__":
    cli()
