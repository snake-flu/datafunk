from datafunk.process_gisaid_sequence_data import *

def run(options):
    process_gisaid_sequence_data(input=options.input, output=options.output, omit_file_list=options.exclude,
                                     exclude_uk=options.exclude_uk, exclude_undated=options.exclude_undated)
