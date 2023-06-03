import os
import time
from copy import deepcopy

import click
import pathops
from cffsubr import subroutinize
from fontTools.fontBuilder import FontBuilder
from fontTools.misc.cliTools import makeOutputFileName
from fontTools.pens.cu2quPen import Cu2QuPen
from fontTools.pens.qu2cuPen import Qu2CuPen
from fontTools.pens.recordingPen import DecomposingRecordingPen
from fontTools.pens.t2CharStringPen import T2CharStringPen
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools.subset import Subsetter
from fontTools.ttLib.scaleUpem import scale_upem
from fontTools.ttLib.ttFont import TTFont, TTLibError, newTable


class TrueTypeToCFFOptions(object):
    def __init__(self):
        self.tolerance: float = 1.0
        self.charstrings_source = "qu2cu"
        self.subroutinize = True
        self.remove_glyphs = True
        self.new_upem = None
        self.recalc_timestamp = False
        self.output_dir = None
        self.overwrite = True


class TrueTypeToCFF(object):
    def __init__(self, font: TTFont):
        self.font = font
        self.options = TrueTypeToCFFOptions()

    def run(self):

        self.font.recalcTimestamp = self.options.recalc_timestamp

        # Prepare the font for conversion
        decomponentize(font=self.font)
        if self.options.new_upem is not None:
            scale_upem(self.font, new_upem=self.options.new_upem)
        if self.options.remove_glyphs:
            remove_unneeded_glyphs(self.font)

        if self.options.charstrings_source == "qu2cu":
            charstrings, glyphs_to_redraw = self.get_qu2cu_charstrings(font=self.font, tolerance=self.options.tolerance)
        else:
            charstrings, glyphs_to_redraw = self.get_t2_charstrings(font=self.font)

        cff_font_info = self.get_cff_font_info()
        post_values = self.get_post_values()

        fb = FontBuilder(font=self.font)
        fb.isTTF = False
        for table in ["glyf", "cvt ", "loca", "fpgm", "prep", "gasp", "LTSH", "hdmx"]:
            if table in fb.font:
                del fb.font[table]
        fb.setupCFF(
            psName=self.font["name"].getDebugName(6),
            charStringsDict=charstrings,
            fontInfo=cff_font_info,
            privateDict={},
        )

        fb.setupDummyDSIG()
        fb.setupMaxp()
        fb.setupPost(**post_values)

        print(glyphs_to_redraw)

        if len(glyphs_to_redraw) > 0:
            temp_source_font: TTFont = deepcopy(self.font)
            from fontTools.subset import Subsetter
            subsetter = Subsetter()
            subsetter.glyph_names_requested = glyphs_to_redraw
            subsetter.options.notdef_glyph = True
            subsetter.subset(temp_source_font)

            generic_error_message(temp_source_font.getGlyphOrder())

            temp_ttf2otf_converter = TrueTypeToCFF(temp_source_font)
            temp_ttf2otf_converter.options.charstrings_source = "t2"
            temp_cff_font = temp_ttf2otf_converter.run()

            temp_otf2ttf_converter = CFFToTrueType(temp_cff_font)
            temp_ttf_font = temp_otf2ttf_converter.run()
            redrawn_charstrings, _ = self.get_qu2cu_charstrings(temp_ttf_font)

            for k, v in redrawn_charstrings.items():
                charstrings[k] = v

            fb.setupCFF(
                psName=self.font["name"].getDebugName(6),
                charStringsDict=charstrings,
                fontInfo=cff_font_info,
                privateDict={},
            )

        if self.options.subroutinize:
            # cffsubr doesn't work with woff/woff2 fonts, we need to temporary set flavor to None
            flavor = fb.font.flavor
            if flavor is not None:
                fb.font.flavor = None
            subroutinize(fb.font)
            fb.font.flavor = flavor

        return fb.font

    def get_cff_font_info(self) -> dict:
        """
        Setup CFF topDict

        :return: A dictionary of the font info.
        """

        font_revision = str(round(self.font["head"].fontRevision, 3)).split(".")
        major_version = str(font_revision[0])
        minor_version = str(font_revision[1]).ljust(3, "0")

        cff_font_info = dict(
            version=".".join([major_version, str(int(minor_version))]),
            FullName=self.font["name"].getBestFullName(),
            FamilyName=self.font["name"].getBestFamilyName(),
            ItalicAngle=self.font["post"].italicAngle,
            UnderlinePosition=self.font["post"].underlinePosition,
            UnderlineThickness=self.font["post"].underlineThickness,
            isFixedPitch=False if self.font["post"].isFixedPitch == 0 else True,
        )

        return cff_font_info

    def get_post_values(self) -> dict:
        post_info = dict(
            italicAngle=round(self.font["post"].italicAngle),
            underlinePosition=self.font["post"].underlinePosition,
            underlineThickness=self.font["post"].underlineThickness,
            isFixedPitch=self.font["post"].isFixedPitch,
            minMemType42=self.font["post"].minMemType42,
            maxMemType42=self.font["post"].maxMemType42,
            minMemType1=self.font["post"].minMemType1,
            maxMemType1=self.font["post"].maxMemType1,
        )
        return post_info

    @staticmethod
    def get_qu2cu_charstrings(font: TTFont, tolerance: float = 1.0) -> (dict, list):
        charstrings = {}
        glyph_set = font.getGlyphSet()
        glyphs_to_redraw = []

        for k, v in glyph_set.items():
            # Correct contours direction and remove overlaps with pathops
            pathops_path = pathops.Path()
            pathops_pen = pathops_path.getPen(glyphSet=glyph_set)
            try:
                glyph_set[k].draw(pathops_pen)
                pathops_path.simplify()
            except TypeError as e:
                generic_warning_message(f"{k}: {e}")
                glyphs_to_redraw.append(k)
                pass

            t2_pen = T2CharStringPen(v.width, glyphSet=glyph_set)
            qu2cu_pen = Qu2CuPen(t2_pen, max_err=tolerance, all_cubic=True, reverse_direction=True)
            try:
                pathops_path.draw(qu2cu_pen)
            except NotImplementedError as e:
                generic_warning_message(f"{k}: {e}")
                qu2cu_pen.all_cubic = False
                pathops_path.draw(qu2cu_pen)
                glyphs_to_redraw.append(k)

            charstring = t2_pen.getCharString()
            charstrings[k] = charstring

        return charstrings, glyphs_to_redraw

    @staticmethod
    def get_t2_charstrings(font: TTFont) -> (dict, list):
        """
        Get CFF charstrings using T2CharStringPen

        :return: CFF charstrings.
        """
        charstrings = {}
        glyph_set = font.getGlyphSet()

        for k, v in glyph_set.items():
            # Draw the glyph with T2CharStringPen and get the charstring
            t2_pen = T2CharStringPen(v.width, glyphSet=glyph_set)
            # qu2cu_pen = Qu2CuPen(t2_pen, max_err=1.0, all_cubic=True, reverse_direction=True)
            glyph_set[k].draw(t2_pen)
            charstring = t2_pen.getCharString()
            charstrings[k] = charstring

        return charstrings, []


class TrueTypeToCFFRunner(object):
    def __init__(self, fonts: list[TTFont]):
        super().__init__()
        self.fonts = fonts
        self.options = TrueTypeToCFFOptions()
        self.is_killed = False

    def run(self) -> None:
        converted_files_count = 0
        start_time = time.time()

        for count, font in enumerate(self.fonts, start=1):
            t = time.time()

            # try:
            print()
            generic_info_message(f"Converting file {count} of {len(self.fonts)}: {font.reader.file.name}")

            # Run the converter
            ttf2otf_converter = TrueTypeToCFF(font=font)
            ttf2otf_converter.options.tolerance = self.options.tolerance / 1000 * font["head"].unitsPerEm
            ttf2otf_converter.options.subroutinize = self.options.subroutinize
            ttf2otf_converter.options.new_upem = self.options.new_upem
            ttf2otf_converter.options.remove_glyphs = self.options.remove_glyphs
            ttf2otf_converter.options.recalc_timestamp = self.options.recalc_timestamp
            ttf2otf_converter.options.charstrings_source = "qu2cu"
            cff_font = ttf2otf_converter.run()

            ext = ".otf" if font.flavor is None else "." + str(font.flavor)
            # Suffix is necessary because without it source woff and woff2 files would be overwritten
            suffix = "" if font.flavor is None else ".otf"
            output_file = makeOutputFileName(
                font.reader.file.name,
                suffix=suffix,
                extension=ext,
                outputDir=self.options.output_dir,
                overWrite=self.options.overwrite,
            )
            cff_font.save(output_file)

            generic_info_message(f"Elapsed time: {round(time.time() - t, 3)} seconds")
            file_saved_message(output_file)
            converted_files_count += 1

            # except Exception as e:
            #     generic_error_message(e)

        print()
        generic_info_message(f"Total files     : {len(self.fonts)}")
        generic_info_message(f"Converted files : {converted_files_count}")
        generic_info_message(f"Elapsed time    : {round(time.time() - start_time, 3)} seconds")


class CFFToTrueTypeOptions(object):
    def __init__(self):
        super().__init__()
        self.max_err = 1.0
        self.post_format = 2.0
        self.reverse_direction = True


class CFFToTrueType(object):
    def __init__(self, font: TTFont):
        self.font = font
        self.options = CFFToTrueTypeOptions()

    def run(self):
        if self.font.sfntVersion != "OTTO":
            raise TTLibError("Not a OpenType font (bad sfntVersion)")
        assert "CFF " in self.font

        glyphOrder = self.font.getGlyphOrder()

        self.font["loca"] = newTable("loca")
        self.font["glyf"] = glyf = newTable("glyf")
        glyf.glyphOrder = glyphOrder
        glyf.glyphs = self.glyphs_to_quadratic(glyphs=self.font.getGlyphSet())
        del self.font["CFF "]
        if "VORG" in self.font:
            del self.font["VORG"]
        glyf.compile(self.font)
        self.update_hmtx(glyf)

        self.font["maxp"] = maxp = newTable("maxp")
        maxp.tableVersion = 0x00010000
        maxp.maxZones = 1
        maxp.maxTwilightPoints = 0
        maxp.maxStorage = 0
        maxp.maxFunctionDefs = 0
        maxp.maxInstructionDefs = 0
        maxp.maxStackElements = 0
        maxp.maxSizeOfInstructions = 0
        maxp.maxComponentElements = max(
            len(g.components if hasattr(g, "components") else []) for g in glyf.glyphs.values()
        )
        maxp.compile(self.font)

        post = self.font["post"]
        post.formatType = self.options.post_format
        post.extraNames = []
        post.mapping = {}
        post.glyphOrder = glyphOrder
        try:
            post.compile(self.font)
        except OverflowError:
            post.formatType = 3

        self.font.sfntVersion = "\000\001\000\000"
        return self.font

    def update_hmtx(self, glyf):
        hmtx = self.font["hmtx"]
        for glyphName, glyph in glyf.glyphs.items():
            if hasattr(glyph, "xMin"):
                hmtx[glyphName] = (hmtx[glyphName][0], glyph.xMin)

    def glyphs_to_quadratic(self, glyphs):
        quadGlyphs = {}
        for gname in glyphs.keys():
            glyph = glyphs[gname]
            ttPen = TTGlyphPen(glyphs)
            cu2quPen = Cu2QuPen(ttPen, max_err=self.options.max_err, reverse_direction=self.options.reverse_direction)
            glyph.draw(cu2quPen)
            quadGlyphs[gname] = ttPen.glyph()
        return quadGlyphs


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
