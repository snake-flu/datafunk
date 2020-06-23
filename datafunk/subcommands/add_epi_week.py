from datafunk.add_epi_week import *

def run(options):
    add_epi_week_column(in_metadata = options.input_metadata,
                 out_metadata = options.output_metadata,
                 date_column = options.date_column,
                 epi_week_column_name = options.epi_week_column_name,
                 epi_day_column_name = options.epi_day_column_name
                )
