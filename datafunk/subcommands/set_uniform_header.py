from datafunk.set_uniform_header import *

def run(options):
    set_uniform_header(options.input_fasta,
                       options.input_metadata,
                       options.output_fasta,
                       options.output_metadata,
                       options.gisaid,
                       options.cog_uk,
                       options.log_file)