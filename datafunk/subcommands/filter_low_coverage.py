from datafunk.filter_low_coverage import *

def run(options):
    remove_low_coverage_sequences(options.input_file, options.output_file, options.threshold)
