from datafunk.add_epi_week import *

def run(options):
    add_epi_week_column(options.input_metadata,
                 options.output_metadata,
                 options.date_column,
                 options.epi_column_name
                )