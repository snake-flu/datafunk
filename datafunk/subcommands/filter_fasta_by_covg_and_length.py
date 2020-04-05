from datafunk.filter_fasta_by_covg_and_length import *

def run(options):
    filter_sequences(options.input_file, options.output_file, options.min_covg, options.min_length)
