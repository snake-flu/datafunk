import pandas as pd
from epiweeks import Week,Year
import datetime


def date_string_to_epi_week(date_string, weeks):
    date_bits = date_string.split("-")
    if len(date_bits) != 3:
        return None
    date = datetime.date(int(date_bits[0]), int(date_bits[1]), int(date_bits[2]))

    for week in weeks:
        if date in week:
            week_string = str(week)
            if "2019" in week_string:
                week_number = "0"
            else:
                week_number = week_string.lstrip("2020")
            return week_number
    return None

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
