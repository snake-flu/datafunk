import pandas as pd
from epiweeks import Week,Year
import datetime
from datetime import datetime
import re
from itertools import chain


def date_string_to_epi_week(date_string):
    """
    parse a date string in YYYY-MM-DD format and return
    cumulative epi week which is cumulative total epidemiological
    weeks since 2019-12-22. Week beginning 2019-12-22 is week 0
    """
    try:
        date = datetime.strptime(date_string, '%Y-%m-%d').date()
    except:
        return ""
    # this is epi-week:
    week = Week.fromdate(date)
    if week.year < 2019 or (week.year == 2019 and week.week < 52):
        return ""
    elif week.year == 2019:
        return("0")
    else:
        cum_epi_week = week.week + len(list(chain(*[[x for x in Year(y).iterweeks()] for y in range(2020, week.year)])))
        return str(cum_epi_week)


def date_string_to_epi_day(date_string):
    """
    parse a date string in YYYY-MM-DD format and return
    cumulative epi day which is cumulative total days since 2019-12-22
    """
    try:
        date = datetime.strptime(date_string, '%Y-%m-%d').date()
    except:
        return ""
    # this is epi-week week:
    week = Week.fromdate(date)
    # this is day 1 of epi-week 0:
    day_one = datetime.strptime("2019-12-22", '%Y-%m-%d').date()
    if week.year < 2019 or (week.year == 2019 and week.week < 52):
        return ""
    else:
        cum_epi_day = (date - day_one).days + 1
        return str(cum_epi_day)


def load_dataframe(metadata_file):
    sep = ','
    if metadata_file.endswith('tsv'):
        sep = '\t'
    df = pd.read_csv(metadata_file, sep=sep)
    return df


def add_epi_week_column(in_metadata, out_metadata, date_column,
                        epi_week_column_name="edin_epi_week",
                        epi_day_column_name=None):

    metadata = load_dataframe(in_metadata)

    epi_week_column = []
    epi_day_column = []
    for i,row in metadata.iterrows():
        date_string = row[date_column]
        epi_week = date_string_to_epi_week(date_string)
        epi_day = date_string_to_epi_day(date_string)
        epi_week_column.append(epi_week)
        epi_day_column.append(epi_day)

    if epi_week_column_name in metadata.columns:
        metadata[epi_week_column_name].update(pd.Series(epi_week_column))
    else:
        metadata[epi_week_column_name] = epi_week_column

    if epi_day_column_name:
        if epi_day_column_name in metadata.columns:
            metadata[epi_day_column_name].update(pd.Series(epi_day_column))
        else:
            metadata[epi_day_column_name] = epi_day_column

    metadata.to_csv(out_metadata, index=False)
