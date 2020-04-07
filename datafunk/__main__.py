import argparse
import sys

import datafunk
import datafunk.subcommands

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

    # _________________________________ repair_names ____________________________#
    subparser_repair_names = subparsers.add_parser("repair_names",
        aliases=['rename_dat_punk'],
        usage="datafunk repair_names --fasta <fasta> --tree <tree> --out <outfile>",
        help="Returns iqtree taxa names to their former glory",
    )

    subparser_repair_names.add_argument("--fasta", action="store", type=str, dest="fasta")
    subparser_repair_names.add_argument("--tree", action="store", type=str, dest="tree")
    subparser_repair_names.add_argument("--out", action="store", type=str, dest="out")

    subparser_repair_names.set_defaults(func=datafunk.subcommands.repair_names.run)

    # ___________________________________________________________________________#

    # _________________________________ remove_fasta ____________________________#
    subparser_remove_fasta = subparsers.add_parser(
        "remove_fasta",
        aliases=['remove_dat_fasta'],
        usage="datafunk remove_fasta -i <input_fasta_file> -f <filter_file> -o <output_file>",
        help="Removing sequences based on a filtering file",
    )

    subparser_remove_fasta.add_argument(
        "-i",
        "--input_file",
        dest="input_file",
        action="store",
        type=str,
        help="Input file: something about the input file format",
    )

    subparser_remove_fasta.add_argument(
        "-f",
        "--filter_file",
        dest="filter_file",
        action="store",
        type=str,
        help="Filter file for filtering based on filter file",
    )

    subparser_remove_fasta.add_argument(
        "-o",
        "--output_file",
        dest="output_file",
        action="store",
        type=str,
        help="Output file name for resulting filtered fasta file",
    )

    subparser_remove_fasta.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )

    subparser_remove_fasta.set_defaults(func=datafunk.subcommands.remove_fasta.run)

    # ___________________________________________________________________________#

    # _________________________________ clean_names ____________________________#
    subparser_clean_names = subparsers.add_parser(
        "clean_names",
        aliases=['clean_dat_name'],
        usage="datafunk clean_names -i <input_metafile> -t <trait> -o <output_file>",
        help="Standardizing country names based on trait given",
    )

    subparser_clean_names.add_argument(
        "-i",
        "--input_metafile",
        dest="input_metafile",
        action="store",
        type=str,
        help="Input file: metafile (csv) for location curation",
    )

    subparser_clean_names.add_argument(
        "-t",
        "--trait",
        dest="trait",
        action="store",
        type=str,
        help="Column name containing the locations",
    )

    subparser_clean_names.add_argument(
        "-o",
        "--output_file",
        dest="output_file",
        action="store",
        type=str,
        help="Output file name for resulting filtered metafile",
    )

    subparser_clean_names.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )

    subparser_clean_names.set_defaults(func=datafunk.subcommands.clean_names.run)

    # ___________________________________________________________________________#

        # _________________________________ merge_fasta ____________________________#
    subparser_merge_fasta = subparsers.add_parser(
        "merge_fasta",
        aliases=['merge_dat_fasta'],
        usage="datafunk merge_fasta -f <folder> -i <input_metafile> -o <output_file>",
        help="Merge fasta file into one single file with removal of duplicates",
    )

    subparser_merge_fasta.add_argument(
        "-f",
        "--folder",
        dest="folder",
        action="store",
        type=str,
        help="Folder containing all fasta files needed to be merged",
    )

    subparser_merge_fasta.add_argument(
        "-i",
        "--input_metafile",
        dest="input_metafile",
        action="store",
        type=str,
        help="Input metafile (csv) for validating sequence information",
    )

    subparser_merge_fasta.add_argument(
        "-o",
        "--output_file",
        dest="output_file",
        action="store",
        type=str,
        help="Output for merged fasta file",
    )

    subparser_merge_fasta.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )

    subparser_merge_fasta.set_defaults(func=datafunk.subcommands.merge_fasta.run)

    # ___________________________________________________________________________#

    # _________________________________ filter_low_coverage ____________________________#
    subparser_filter_low_coverage = subparsers.add_parser(
        "filter_low_coverage",
        aliases=['filter_dat_fasta'],
        usage="datafunk filter_low_coverage -i <input_fasta> -t <threshold> [-o <output_fasta>]",
        help="Removes sequences where the fraction of non-N bases falls below the threshold",
    )

    subparser_filter_low_coverage.add_argument(
        "-i",
        "--input_file",
        dest="input_file",
        action="store",
        required=True,
        type=str,
        help="Input FASTA",
    )

    subparser_filter_low_coverage.add_argument(
        "-t",
        "--threshold",
        dest="threshold",
        action="store",
        required=True,
        type=int,
        help="Integer representing the percentage threshold. Sequences with coverage (strictly) "
             "less than this will be excluded from the filtered file.",
    )

    subparser_filter_low_coverage.add_argument(
        "-o",
        "--output_file",
        dest="output_file",
        action="store",
        default=None,
        type=str,
        help="Output file name for resulting filtered FASTA (default adds .filtered to input file name)",
    )

    subparser_filter_low_coverage.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )

    subparser_filter_low_coverage.set_defaults(func=datafunk.subcommands.filter_low_coverage.run)

    # ___________________________________________________________________________#

    # _________________________________ process_gisaid_sequence_data ____________________________#
    subparser_process_gisaid_sequence_data = subparsers.add_parser(
        "process_gisaid_sequence_data",
        aliases=['get_new_horrors'],
        usage="datafunk process_gisaid_sequence_data -i <input.json OR input.fasta> [-o <output.fasta>] [-e file1 -e file2 ...] [--stdout]",
        help="Process raw sequence data in fasta or json format",
        description="Process raw sequence data in fasta or json format",
    )

    subparser_process_gisaid_sequence_data._action_groups.pop()
    required_process_gisaid_sequence_data = subparser_process_gisaid_sequence_data.add_argument_group('required arguments')
    optional_process_gisaid_sequence_data = subparser_process_gisaid_sequence_data.add_argument_group('optional arguments')

    required_process_gisaid_sequence_data.add_argument(
        '-i',
        '--input',
        required=True,
        metavar='GISAID.fasta OR GISAID.json',
        help='Sequence data in FASTA/json format'
    )

    optional_process_gisaid_sequence_data.add_argument(
        '-o',
        '--output',
        required=False,
        metavar = 'OUTPUT.fasta',
        help='FASTA format file to write'
    )

    optional_process_gisaid_sequence_data.add_argument(
        '-e',
        '--exclude',
        action='append',
        required=False,
        metavar = 'FILE',
        help='A file that contains (anywhere) EPI_ISL_###### IDs to exclude (can provide many files, '
             'e.g. -e FILE1 -e FILE2 ...)'
    )


    optional_process_gisaid_sequence_data.add_argument(
        '--stdout',
        required=False,
        action='store_true',
        help='Overides -o/--output if present and prints output to stdout',
    )

    subparser_process_gisaid_sequence_data.set_defaults(func=datafunk.subcommands.process_gisaid_sequence_data.run)

    # _________________________________ sam_2_fasta _____________________________#

    subparser_sam_2_fasta = subparsers.add_parser(
        "sam_2_fasta",
        usage="datafunk sam_2_fasta -s <input.sam> -r <reference.fasta> [-o <output.fasta>] [-t [INT]:[INT]] [--prefix_ref] [--stdout]",
        help="Convert sam format alignment to fasta format multiple alignment, with optional trimming",
        description="aligned sam -> fasta (with optional trim to user-defined (reference) co-ordinates)",
    )

    subparser_sam_2_fasta._action_groups.pop()
    required_sam_2_fasta = subparser_sam_2_fasta.add_argument_group('required arguments')
    optional_sam_2_fasta = subparser_sam_2_fasta.add_argument_group('optional arguments')

    required_sam_2_fasta.add_argument(
        '-s', '--sam',
        help='samfile',
        required=True,
        metavar='in.sam'
                        )
    required_sam_2_fasta.add_argument(
        '-r', '--reference',
        help='reference',
        required=True,
        metavar='reference.fasta'
                        )
    optional_sam_2_fasta.add_argument(
        '-o', '--output',
        help='FASTA format file to write.',
        required=False,
        metavar='out.fasta'
                        )
    optional_sam_2_fasta.add_argument(
        '-t', '--trim',
        help='trim the alignment to these coordinates (0-based, half-open)',
        required=False,
        metavar='[[start]:[end]]'
                        )
    optional_sam_2_fasta.add_argument(
        '--prefix_ref',
        help='write the reference sequence at the beginning of the file',
        required=False,
        action='store_true'
                        )
    optional_sam_2_fasta.add_argument(
        '--stdout',
        help='Overides -o/--output if present and prints output to stdout',
        required=False,
        action='store_true'
                        )

    subparser_sam_2_fasta.set_defaults(func=datafunk.subcommands.sam_2_fasta.run)

    # _________________________________ gisaid_json_2_metadata _____________________________#

    subparser_gisaid_json_2_metadata = subparsers.add_parser(
        "gisaid_json_2_metadata",
        usage="json_2_csv.py [-h] -n gisaid.json -c <metadata.csv / False> -o metadata.csv -e omissions.txt",
        help = "Add the info from a Gisaid json dump to an existing metadata table (or create a new one)",
        description="""Add the info from a Gisaid json dump to an existing metadata table (or create a new one)""")

    subparser_gisaid_json_2_metadata._action_groups.pop()
    required_gisaid_json_2_metadata = subparser_gisaid_json_2_metadata.add_argument_group('required arguments')
    optional_gisaid_json_2_metadata = subparser_gisaid_json_2_metadata.add_argument_group('optional arguments')

    required_gisaid_json_2_metadata.add_argument('-n', '--new',
                        help='Most recent Gisaid json dump',
                        required=True,
                        metavar = 'gisaid.json')
    required_gisaid_json_2_metadata.add_argument('-c', '--csv',
                        help='Last metadata table (csv format), or \'False\' if you really want '
                             '(but you will lose date stamp information from previous dumps)',
                        required=True,
                        metavar = '<metadata.csv / False>')
    required_gisaid_json_2_metadata.add_argument('-o', '--output',
                        help='New csv file to write',
                        required=True,
                        metavar='metadata.csv')
    required_gisaid_json_2_metadata.add_argument('-e', '--exclude',
                        action='append',
                        required=True,
                        metavar = 'omissions.txt',
                        help='A file that contains (anywhere) EPI_ISL_###### IDs to exclude (can provide more than one file, '
                             'e.g. -e FILE1 -e FILE2 ...)')


    subparser_gisaid_json_2_metadata.set_defaults(func=datafunk.subcommands.gisaid_json_2_metadata.run)

    # ___________________________________________________________________________#

    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
