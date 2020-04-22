from datafunk.gisaid_json_2_metadata import *

def run(options):
    gisaid_json_2_metadata(json = options.new, \
                           output = options.output-metadata, \
                           args_csv = options.csv, \
                           args_omit_file_list = options.exclude,
                           args_lineages = options.lineages)
