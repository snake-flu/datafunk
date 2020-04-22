from datafunk.process_gisaid_data import *

def run(options):
    process_gisaid_data(input_json=options.json,
                        input_omit_file_list=options.exclude,
                        input_metadata=options.input_metadata,
                        output_fasta=options.output_sequences,
                        output_metadata = options.output_metadata,
                        exclude_uk=options.exclude_uk,
                        exclude_undated=options.exclude_undated,
                        exclude_subsampled = not(options.include_subsampled),
                        exclude_omitted_file = not(options.include_omitted_file))
