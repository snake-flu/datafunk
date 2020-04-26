import pandas as pd
from epiweeks import Week,Year
import datetime
from datetime import datetime
import re


def date_string_to_epi_week(date_string, weeks):
    # check the date:
    regex = re.compile('\d{4}-\d{2}-\d{2}')
    match = re.search(regex, date_string)
    if not match:
        return(None)

    date = datetime.strptime(date_string, '%Y-%m-%d').date()

    week = Week.fromdate(date)

    if week in weeks:
        if str(week)[0:4] == '2019':
            return('0')
        else:
            return(str(week.weektuple()[1]))
    else:
        return(None)

def load_dataframe(metadata_file):
    sep = ','
    if metadata_file.endswith('tsv'):
        sep = '\t'
    df = pd.read_csv(metadata_file, sep=sep)
    return df

def add_epi_week_column(in_metadata, out_metadata, date_column, epi_column_name="edin_epi_week"):
    metadata = load_dataframe(in_metadata)

    last_2019 = Week(2019, 52)
    weeks = list(Year(2020).iterweeks())
    weeks.append(last_2019)

    epi_column = []
    for i,row in metadata.iterrows():
        date_string = row[date_column]
        epi_week = date_string_to_epi_week(date_string, weeks)
        epi_column.append(epi_week)
    metadata[epi_column_name] = epi_column

    metadata.to_csv(out_metadata, index=False)
