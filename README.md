# TTF2OTF

Converts TTF fonts (or TrueType flavored woff/woff2 web fonts) to OTF fonts (or CFF flavored woff/woff2 web fonts).

**Usage**:

    ttf2otf.py [OPTIONS] INPUT_PATH

**Options**:

      -t, --tolerance FLOAT RANGE   Conversion tolerance (0-3, default 1). Low
                                    tolerance adds more points but keeps shapes.
                                    High tolerance adds few points but may change
                                    shape.  [0<=x<=3]
      --scale-upm INTEGER RANGE     Scale units-per-em of source font to passed
                                    value before converting to OTF  [1000<=x<8192]
      --remove-glyphs               Remove NULL and CR glyphs from output fonts
      --no-subr                     Do not subroutinize converted fonts with
                                    cffsubr
      -out, --output-dir DIRECTORY  Specify the directory where output files are
                                    to be saved. If output_dir doesn't exist, will
                                    be created. If not specified, files are saved
                                    to the source folder.
      --recalc-timestamp            Keep the original font 'modified' timestamp
                                    (head.modified) or set it to current time. By
                                    default, original timestamp is kept.
      --no-overwrite                Overwrite existing output files or save them
                                    to a new file (numbers are appended at the end
                                    of file name). By default, files are
                                    overwritten.
      --help                        Show this message and exit.
