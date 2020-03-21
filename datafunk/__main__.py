import argparse
import sys

import datafunk

def main(args=None):
    parser = argparse.ArgumentParser(
        prog="datafunk",
        usage="datafunk <subcommand> <options>",
        description="Miscellaneous data manipulation tools",
    )

    parser.add_argument("--version", action="version", version=datafunk.__version__)
    subparsers = parser.add_subparsers(
        title="Available subcommands", help="", metavar=""
    )

    #_____________________________ remove_dat_junk ______________________________#
    subparser_remove_dat_junk = subparsers.add_parser(
        "remove_dat_junk",
        usage="datafunk remove_dat_junk -i <input>",
        help="Example command",
    )

    subparser_remove_dat_junk.add_argument(
        "-i",
        "--input_file",
        dest="input_file",
        action="store",
        type=str,
        help="Input file: something about the input file format",
    )

    subparser_remove_dat_junk.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )

    subparser_remove_dat_junk.set_defaults(func=datafunk.subcommands.remove_dat_junk.run)

    # ___________________________________________________________________________#

    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
