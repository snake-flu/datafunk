import pandas as pd
import numpy as np
import re
from Bio import SeqIO

def strip_nasties(name):
    return name.replace(" ","_")\
        .replace("hCoV-19/","")\
        .replace("hCov-19/","")\
        .replace("PENDING", "")\
        .replace("UPLOADED", "")\
        .replace("None", "")\
        .replace("UK-ENG", "England")\
        .replace("UK-SCT", "Scotland")\
        .replace("UK-WLS", "Wales")

def load_dataframe(metadata_file):
    sep = ','
    if metadata_file.endswith('tsv'):
        sep = '\t'
    df = pd.read_csv(metadata_file, sep=sep)
    return df

def add_header_column(df, columns):
    headers = []
    for i,row in df.iterrows():
        fields = [row[c] if isinstance(row[c], str) else "" for c in columns ]
        header = '|'.join(fields)
        header = strip_nasties(header)
        headers.append(header)
    df['header'] = headers
    return df

def parse_virus_name(header):
    regex = r"[A-Za-z _]+/[\w_-]+/\d+"
    regex = re.compile(regex)
    match = re.search(regex, header)
    if match:
        return match.group()
    else:
        return ""

def parse_country(header):
    regex = r"\|[A-Za-z _]+[\|\/]"
    regex = re.compile(regex)
    match = re.search(regex, header)
    if match:
        return match.group()[1:-1]
    else:
        return ""

def parse_country_from_virus_name(name):
    return name.split('/')[0]

def parse_date(header):
    regex = r"\|\d+-\d+-\d+"
    regex = re.compile(regex)
    match = re.search(regex, header)
    if match:
        return match.group()[1:]
    else:
        return ""

def get_new_header(header):
    name = parse_virus_name(header)
    country = parse_country(header)
    if country == "":
        country = parse_country_from_virus_name(name)
    new_header = name + '|' + country + '|' + parse_date(header)
    new_header = strip_nasties(new_header)
    return new_header

def get_new_header_second_attempt(header):
    name = parse_virus_name(header)
    country = parse_country_from_virus_name(name)
    new_header = name + '|' + country + '|' + parse_date(header)
    new_header = strip_nasties(new_header)
    return new_header

def header_found_in_column(header, df, column):
    return (df[column] == header).any()

def id_found_in_column(header, df, column):
    sample_id = header.split('|')[0].split('/')[1]
    if sample_id == "":
        return False
    return header_found_in_column(sample_id, df, column)

def update_df_if_id_found_in_column(header, df, column):
    sample_id = header.split('|')[0].split('/')[1]
    if sample_id == "":
        return df
    if (df[column] == sample_id).any():
        df.loc[df[column] == sample_id,"header"] = header
    return df

def header_duplicated_in_column(header, df, column):
    return df[df[column] == header][column].duplicated().any()

#GISAID: covv_virus_name, edin_admin_0 or covv_location.split(' / ')[1], covv_collection_date
gisaid_columns = ['covv_virus_name', 'edin_admin_0', 'covv_collection_date']

#COG-UK: secondary_accession adm1 collection_date
coguk_columns = ['secondary_identifier', 'adm1', 'collection_date']

def set_uniform_header(input_fasta, input_metadata, output_fasta, output_metadata, gisaid, cog_uk, log_file):

    metadata = load_dataframe(input_metadata)
    if gisaid:
        metadata = add_header_column(metadata, gisaid_columns)
    elif cog_uk:
        metadata = add_header_column(metadata, coguk_columns)
    else:
        sys.exit("Must use either --gisaid or --cog_uk flag")

    if log_file:
        log_handle = open(log_file, "w")
    else:
        log_handle = sys.stdout

    with open(input_fasta) as in_fasta, open(output_fasta, 'w') as out_fasta:
        for record in SeqIO.parse(in_fasta, "fasta"):
            header = get_new_header(record.description)
            if not header_found_in_column(header, metadata, "header"):
                header = get_new_header_second_attempt(record.description)
            if cog_uk and not header_found_in_column(header, metadata, "header") \
                    and id_found_in_column(header, metadata, "central_sample_id"):
                metadata = update_df_if_id_found_in_column(header, metadata, "central_sample_id")

            if not header_found_in_column(header, metadata, "header"):
                log_handle.write("Could not find header %s parsed from record %s in metadata table\n" %(header, record.id))
            else:
                if header_duplicated_in_column(header, metadata, "header"):
                    log_handle.write("Header %s parsed from record %s has duplicate entries in metadata table\n" % (
                    header, record.id))
                record.id = header
                record.description = ""
                SeqIO.write(record, out_fasta, "fasta-2line")

    metadata.to_csv(output_metadata)

    if log_handle:
        log_handle.close()