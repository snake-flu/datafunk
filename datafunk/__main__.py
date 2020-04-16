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
        default='removed.fasta',
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
        default='cleaned.csv',
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
        default='filtered.fasta',
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

    # _________________________________ filter_fasta_by_covg_and_length ____________________________#
    subparser_filter_fasta_by_covg_and_length = subparsers.add_parser(
        "filter_fasta_by_covg_and_length",
        aliases=['filter_dat_fasta'],
        usage="datafunk filter_fasta_by_covg_and_length -i <input_fasta> -t <threshold> [-o <output_fasta>]",
        help="Removes sequences where the fraction of non-N bases falls below the threshold",
    )

    subparser_filter_fasta_by_covg_and_length.add_argument(
        "-i",
        "--input_file",
        dest="input_file",
        action="store",
        required=True,
        type=str,
        help="Input FASTA",
    )

    subparser_filter_fasta_by_covg_and_length.add_argument(
        "--min_covg",
        dest="min_covg",
        action="store",
        required=False,
        type=int,
        help="Integer representing the minimum coverage percentage threshold. Sequences with coverage "
             "(strictly) less than this will be excluded from the filtered file.",
    )

    subparser_filter_fasta_by_covg_and_length.add_argument(
        "--min_length",
        dest="min_length",
        action="store",
        required=False,
        type=int,
        help="Integer representing the minimum length threshold. Sequences with length (strictly) "
             "less than this will be excluded from the filtered file. Default: ?",
    )

    subparser_filter_fasta_by_covg_and_length.add_argument(
        "-o",
        "--output_file",
        dest="output_file",
        action="store",
        default=None,
        type=str,
        help="Output file name for resulting filtered FASTA (default adds .filtered to input file name)",
    )

    subparser_filter_fasta_by_covg_and_length.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )

    subparser_filter_fasta_by_covg_and_length.set_defaults(func=datafunk.subcommands.filter_fasta_by_covg_and_length.run)

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
        help='FASTA format file to write, print to stdout if unspecified'
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
        '--exclude_uk',
        required=False,
        action='store_true',
        help='Removes all GISAID entries with containing England, Ireland, Scotland or Wales',
    )

    optional_process_gisaid_sequence_data.add_argument(
        '--exclude_undated',
        required=False,
        action='store_true',
        help='Removes all GISAID entries with an incomplete date',
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
        help='FASTA format file to write. Prints to stdout if not specified',
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
        '--log_inserts',
        help='log insertions relative to the reference',
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


    # ___________________________________________________________________________#

    # _________________________________ phylotype_consensus ____________________________#
    subparser_phylotype_consensus = subparsers.add_parser(
        "phylotype_consensus",
        aliases=['phylotype_consensus'],
        usage="datafunk phylotype_consensus -i <input_fasta> -m <input_metafile> -c <clade_file> -o <output_folder>",
        help="Split sequences according to lineage and create an consensus",
    )

    subparser_phylotype_consensus.add_argument(
        "-i",
        "--input_fasta",
        dest="input_fasta",
        action="store",
        type=str,
        help="Fasta file for splitting into phylotypes",
    )

    subparser_phylotype_consensus.add_argument(
        "-m",
        "--input_metafile",
        dest="input_metafile",
        action="store",
        type=str,
        help="Input metafile (csv) with phylotype information",
    )


    subparser_phylotype_consensus.add_argument(
        "-c",
        "--clade_file",
        dest="clade_file",
        action="store",
        type=str,
        help="Clade file stating the phylotypes needed to be grouped",
    )

    subparser_phylotype_consensus.add_argument(
        "-o",
        "--output_folder",
        dest="output_folder",
        action="store",
        default="./",
        type=str,
        help="Output folder for the phylotype fasta files and consensus file",
    )

    subparser_phylotype_consensus.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )

    subparser_phylotype_consensus.set_defaults(func=datafunk.subcommands.phylotype_consensus.run)

    # _________________________________ gisaid_json_2_metadata ____________________________#

    subparser_gisaid_json_2_metadata = subparsers.add_parser(
        """gisaid_json_2_metadata""",
        usage="""datafunk gisaid_json_2_metadata [-h] -n gisaid.json -c <OLD_metadata.csv / False> -o NEW_metadata.csv -e omissions.txt""",
        description="""Add the info from a Gisaid json dump to an existing metadata table (or create a new one)""",
        help="""Add the info from a Gisaid json dump to an existing metadata table (or create a new one)""")

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
                        metavar = '<OLD_metadata.csv / False>')
    required_gisaid_json_2_metadata.add_argument('-o', '--output',
                        help='New csv file to write',
                        required=True,
                        metavar='NEW_metadata.csv')
    optional_gisaid_json_2_metadata.add_argument(
        '-l',
        '--lineages',
        required=False,
        metavar = 'lineages.csv',
        help='csv file of ineages to include'
    )
    required_gisaid_json_2_metadata.add_argument('-e', '--exclude',
                        action='append',
                        required=True,
                        metavar = 'omissions.txt',
                        help='A file that contains (anywhere) EPI_ISL_###### IDs to exclude (can provide more than one file, '
                             'e.g. -e FILE1 -e FILE2 ...)')

    subparser_gisaid_json_2_metadata.set_defaults(func=datafunk.subcommands.gisaid_json_2_metadata.run)

    # _________________________________ set_uniform_header ____________________________#

    subparser_set_uniform_header = subparsers.add_parser(
        "set_uniform_header",
        usage="datafunk set_uniform_header -i <input_fasta> -t <threshold> [-o <output_fasta>]",
        help="Adds a header column to metadata table and converts fasta to have a uniform header",
    )

    subparser_set_uniform_header.add_argument(
        "--input_fasta",
        dest="input_fasta",
        required=True,
        type=str,
        help="Input FASTA",
    )
    subparser_set_uniform_header.add_argument(
        "--input_metadata",
        dest="input_metadata",
        required=True,
        type=str,
        help="Input CSV or TSV",
    )
    subparser_set_uniform_header.add_argument(
        "--output_fasta",
        dest="output_fasta",
        required=True,
        type=str,
        help="Input FASTA",
    )
    subparser_set_uniform_header.add_argument(
        "--output_metadata",
        dest="output_metadata",
        required=True,
        type=str,
        help="Input CSV or TSV",
    )
    subparser_set_uniform_header.add_argument(
        "--gisaid",
        dest="gisaid",
        action="store_true",
        required=False,
        help="Input data is from GISAID",
    )
    subparser_set_uniform_header.add_argument(
        "--cog_uk",
        dest="cog_uk",
        action="store_true",
        required=False,
        help="Input data is from COG-UK",
    )
    subparser_set_uniform_header.add_argument(
        "--log",
        dest="log_file",
        required=False,
        help="Log file to use (otherwise uses stdout)"
    )
    subparser_set_uniform_header.add_argument(
        "--column_name",
        dest="column_name",
        required=False,
        default='sequence_name',
        help="Name of column in metadata corresponding to fasta header"
    )
    subparser_set_uniform_header.add_argument(
        "--index_column",
        dest="index_column",
        required=False,
        help="Name of column in metadata to parse for string matching with fasta header"
    )

    subparser_set_uniform_header.set_defaults(func=datafunk.subcommands.set_uniform_header.run)

    # _________________________________ add_epi_week ____________________________#

    subparser_add_epi_week = subparsers.add_parser(
        "add_epi_week",
        usage="datafunk add_epi_week -i <input_metadata> -o <output_metadata> --date_column <column> [--epi_column_name <column>]",
        help="Adds a column for epi week to metadata table",
    )

    subparser_add_epi_week.add_argument(
        "-i", "--input_metadata",
        dest="input_metadata",
        required=True,
        type=str,
        help="Input CSV or TSV",
    )
    subparser_add_epi_week.add_argument(
        "-o", "--output_metadata",
        dest="output_metadata",
        required=True,
        type=str,
        help="Input CSV or TSV",
    )
    subparser_add_epi_week.add_argument(
        "--date_column",
        dest="date_column",
        required=True,
        help="Column name to convert to epi week",
    )
    subparser_add_epi_week.add_argument(
        "--epi_column_name",
        dest="epi_column_name",
        required=False,
        help="Column name for epi week column",
    )

    subparser_add_epi_week.set_defaults(func=datafunk.subcommands.add_epi_week.run)

    # ___________________________________________________________________________#

    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
