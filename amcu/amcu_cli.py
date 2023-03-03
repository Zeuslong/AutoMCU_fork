import os
import argparse
from automcu import AutoMCU


def main():
    """
    Required and optional arguments to be included in the CLI, while running the code.
    The arguments are defined in the main function. The arguments are parsed in the __main__.py file, or used in the run function in the automcu.py file.
    """
    parser = argparse.ArgumentParser(description="Run AutoMCU on an image")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Force overwrite of existing output file",
    )
    parser.add_argument(
        "--config",
        "-c",
        default=None,
        help="Get settings from JSON-formatted config file,"
        " instead of needing to specify command line args",
    )
    parser.add_argument(
        "--input_hdr",
        default=None,
        help="Specify name of input ENVI hdr in case it is not"
        " the standard path.hdr or splitext(path).hdr",
    )
    parser.add_argument(
        "--wl_range",
        "-w",
        default=[],
        action="append",
        help="A 2-element list or tuple of floats, i.e."
        " (600.0,750), that will be used to identify"
        " the band ranges used for unmixing."
        " Data from the bands specified by"
        " these wavelengths will be lumped into a single"
        " array for fitting. A minimum of one wl window"
        " must be specified. Both wl_range and band_range"
        " cannot be specified",
    )
    parser.add_argument(
        "--band_range",
        "-b",
        default=[],
        action="append",
        help="A 2-element list or tuple of integers, i.e."
        " (65,72), that will be used to identify"
        " the band ranges used for unmixing."
        " Band numbers are 0-based, so first band is 0."
        " Data from the bands specified by"
        " these windows will be lumped into a single"
        " array for fitting. A minimum of one band window"
        " must be specified. Both wl_range and band_range"
        " cannot be specified",
    )
    parser.add_argument(
        "--notie",
        dest="tied",
        action="store_false",
        help="Reflectance spectra from each band range are"
        " subtracted from the first window value to"
        " help relieve brightness issues. (First value"
        " in unmixing array becomes 0.0 for each"
        " window's start). By default this is on. Use this"
        " flag to disable.",
    )
    parser.add_argument(
        "--divide",
        action="store_true",
        help="Reflectance spectra from each band range are"
        " divided by the first window value to"
        " help relieve brightness issues. (First value"
        " in unmixing array becomes 1.0 for each"
        " window's start). By default this is off.",
    )
    parser.add_argument(
        "--shade", action="store_true", help="Add an endmember for shade"
    )
    parser.add_argument(
        "--sum_to_one",
        "-1",
        action="store_true",
        help="Add an element with value 1.0 to the end of the"
        " unmixing array and the input refl to help"
        " coerce sums of the endmember fractions to be"
        " near 1.0.",
    )
    parser.add_argument(
        "--iterations",
        "--count",
        type=int,
        default=30,
        help="How many MC runs to do for each pixel. " " Default 30",
    )

    parser.add_argument(
        "--num_blocks",
        default=[],
        action="append",
        help="How many blocks to unmix. "
        " Default all"
        "A 2-element list or tuple of int, i.e."
        " (100,200), that will be used to identify"
        " the blocks ranges used for unmixing.",
    )
    ##define a argument having a range for a parameter --num_blocks

    parser.add_argument(
        "--outscale",
        type=float,
        default=1000,
        help="Scale value applied to unmixing results to make"
        " them fit into data type specified for output.",
    )
    parser.add_argument(
        "--nodata",
        type=float,
        default=None,
        help="Output no data value written into metadata. If"
        " default None, missing data in output will be 0,"
        " but this will not be recorded in metadata.",
    )
    parser.add_argument(
        "--scale",
        type=float,
        default=1.0,
        help="Apply this scale to input data to get" " reflectance scaled from 0 to 1",
    )
    parser.add_argument(
        "--dtype",
        "--of",
        default="int16",
        help="Rasterio-usable specification of output map" " data type.",
    )
    parser.add_argument(
        "--otype",
        "--ot",
        default="GTiff",
        help="GDAL shortname for output data format.",
    )
    parser.add_argument(
        "--co",
        "-o",
        default=[],
        action="append",
        help="GDAL data format writing options.",
    )
    parser.add_argument(
        "--emfwhm",
        action="store_true",
        help="Assume second column of em csvs are fwhm values."
        " By default, average inter-band spacing is used",
    )
    parser.add_argument(
        "--nointerp",
        action="store_true",
        help="Assume em csvs have same count and wavelengths"
        " as input reflectance data (will be checked)",
    )
    parser.add_argument(
        "--match_c",
        action="store_true",
        help="Recreate a bug in original C code that shifted"
        " input reflectance data by one index/band",
    )
    parser.add_argument(
        "--names",
        default=[],
        action="append",
        metavar="STR[,STR]",
        help="Supply names for the classes defined by the"
        " given CSV files. Can be specified multiple"
        " times or using comma separated values",
    )
    parser.add_argument("input", help="ENVI-formatted reflectance image")
    parser.add_argument("output", help="AutoMCU result map filename")
    parser.add_argument(
        "emlist",
        nargs="+",
        metavar="CSV<:INT>",
        help="CSV file with list of samples for a given"
        " endmember class. CSV should should be arranged"
        " such that samples are columns, and bands are"
        " rows, with the first column containing the band"
        " center wavelengths in nanometers. At least one"
        " class must be specified. If a : and integer n"
        " immediately follow the filename, then only the"
        " first n samples will be read",
    )
    args = parser.parse_args()
    if args.verbose:
        print("Args:")
        for k, v in vars(args).items():
            print(f"{k:12s}: {v}")

    ##Check that we have ems
    em_csvs = []
    em_counts = []
    for st in args.emlist:
        parts = st.strip().split(":")
        if len(parts) > 1:
            em_csvs.append(parts[0])
            if len(parts[1]) > 0:
                em_counts.append(int(parts[1]))
            else:
                em_counts.append(0)
        else:
            em_csvs.append(parts[0])
            em_counts.append(0)

    em_names = []
    if args.names:
        for n in args.names:
            parts = n.strip().split(",")
            em_names.extend(parts)

    ###########################
    ##Process band or wl ranges
    ###########################
    ##Check that we have band or wavelength windows
    if (len(args.band_range) > 0) and (len(args.wl_range) > 0):
        raise RuntimeError("Can't mix wl_range and band_range entries")

    ##Convert to floats or ints
    if args.band_range:
        band_ranges = []
        wl_ranges = None
        for bstr in args.band_range:
            nobrackets = bstr.strip("[](){}")
            parts = nobrackets.split(",", maxsplit=2)
            if len(parts) < 2:
                raise RuntimeError(
                    "Could not convert band range to integers:" f" {bstr}"
                )
            band_ranges.append([int(p) for p in parts])

    elif args.wl_range:
        band_ranges = None
        wl_ranges = []
        for wstr in args.wl_range:
            nobrackets = wstr.strip("[](){}")
            parts = nobrackets.split(",", maxsplit=2)
            if len(parts) < 2:
                raise RuntimeError("Could not convert wl range to integers:" f" {wstr}")
            wl_ranges.append([float(p) for p in parts])
    else:
        raise RuntimeError(
            "Band or wavelength ranges must be specified with"
            "either --wl_range or --band_range"
        )

    if args.num_blocks:
        nblocks = []
        for nstr in args.num_blocks:
            nobrackets = nstr.strip("[](){}")
            parts = nobrackets.split(",", maxsplit=2)
            if len(parts) < 2:
                raise RuntimeError(
                    "Could not convert number of blocks to integers:" f" {nstr}"
                )
            nblocks.append([int(p) for p in parts])
    ##As a failsafe to prevent writing over one of the input files if command
    ## line argument order is mixed up, don't write output if the file exists
    if os.path.exists(args.output) and (not args.overwrite):
        raise RuntimeError(
            f"File {args.output} already exists. Check args"
            "with -v or use --overwrite if this is intentional"
        )

    #############################
    ##Create the AutoMCU instance
    #############################
    amcu = AutoMCU(
        em_csvs,
        em_counts,
        em_names=em_names,
        band_ranges=band_ranges,
        wl_ranges=wl_ranges,
        tied=args.tied,
        divide=args.divide,
        shade=args.shade,
        sum_to_one=args.sum_to_one,
        iterations=args.iterations,
        outscale=args.outscale,
        nodata=args.nodata,
        dtype=args.dtype,
        otype=args.otype,
        co=args.co,
        emfwhm=args.emfwhm,
        nointerp=args.nointerp,
        match_c=args.match_c,
    )

    ##########################
    ##Run AutoMCU on the image
    ##########################
    amcu.apply_image(
        args.input,
        args.scale,
        args.output,
        hdrpath=args.input_hdr,
        nblocks=args.num_blocks,
    )

    return


if __name__ == "__main__":
    main()
