from datafunk.add_header_column import *

def run(options):
    add_header_column(options.input_fasta,
                       options.input_metadata,
                       options.output_fasta,
                       options.output_metadata,
                       options.log_file,
                       options.column_name,
                       options.columns,
                       options.gisaid,
                       options.cog_uk)

