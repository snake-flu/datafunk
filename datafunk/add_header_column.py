import pandas as pd
import numpy as np
import re
from Bio import SeqIO
import sys

def load_dataframe(metadata_file):
    sep = ','
    if metadata_file.endswith('tsv'):
        sep = '\t'
    na_values = ["None", ""]
    df = pd.read_csv(metadata_file, sep=sep, na_values=na_values)
    return df

def parse_virus_name(header):
    regex = r"[A-Za-z _-]+/[\w_-]+/[A-Z:_]*\d+"
    regex = re.compile(regex)
    match = re.search(regex, header)
    if match:
        return match.group().split('/')[1]
    else:
        return None

def add_header_entry(df, header_name, columns, log_handle, column_name='header', count=0):
    sample_id = parse_virus_name(header_name)
    if sample_id is None:
        log_handle.write("Could not parse id from %s\n" % header_name)
        return

    indexes = []
    for column in columns:
        if (df[column] == sample_id).any():
            df_subset = df[df[column] == sample_id]
            for i,row in df_subset.iterrows():
                found = True
                for column in columns:
                    if "/"+str(row[column])+"/" not in header_name \
                        and "|"+str(row[column]) not in header_name:
                        found = False
                        break
                if found:
                    indexes.append(i)
    if len(indexes) == 0 and 'secondary_identifier' in df.columns:
        for i, row in enumerate(df['secondary_identifier'].values):
            if sample_id in str(row):
                indexes.append(i)
        if len(indexes) > 1:
            df_subset = df.loc[indexes,:]
            indexes = []
            for i, row in df_subset.iterrows():
                found = True
                for column in columns[1:]:
                    if "/" + str(row[column]) + "/" not in header_name \
                            and "|" + str(row[column]) not in header_name:
                        found = False
                        break
                if found:
                    indexes.append(i)

    if len(indexes) > 10:
        log_handle.write(
            "Found %i rows corresponding to sample_id %s and header %s which must be wrong\n" % (len(indexes),
                                                                                            sample_id, header_name))
        indexes = []
    elif len(indexes) == 0:
        log_handle.write("Could not find corresponding row for %s\n" % header_name)
        return
    if count > 0 and len(indexes) > count:
        indexes = indexes[count:]
        header_name += "_" + str(count)
    for index in indexes:
        if df.loc[index, column_name] != None:
            log_handle.write("Overwriting %s with %s\n" % (df.loc[index, column_name], header_name))
    df.loc[indexes, column_name] = header_name
    return

#GISAID: covv_virus_name, edin_admin_0 or covv_location.split(' / ')[1], covv_collection_date
gisaid_columns = ['covv_virus_name', 'edin_admin_0', 'covv_collection_date', 'covv_accession_id']

#COG-UK: secondary_accession adm1 collection_date
coguk_columns = ['central_sample_id', 'sequencing_org_code', "sequencing_submission_date"]

def add_header_column(input_fasta, input_metadata, output_fasta, output_metadata, log_file, column_name, columns,
                      gisaid, cog_uk):

    metadata = load_dataframe(input_metadata)
    metadata[column_name] = None

    if log_file:
        log_handle = open(log_file, "w")
    else:
        log_handle = sys.stdout

    with open(input_fasta) as in_fasta, open(output_fasta, 'w') as out_fasta:
        found_headers = []
        for record in SeqIO.parse(in_fasta, "fasta"):
            header_name = record.description
            count = 0
            if len(header_name) == 0:
                log_handle.write("Bad header %s in input fasta\n" % (record.id))
                continue
            elif header_name in found_headers:
                count = found_headers.count(header_name)
            found_headers.append(header_name)

            if gisaid:
                add_header_entry(metadata, header_name, gisaid_columns, log_handle, column_name, count)
            elif cog_uk:
                add_header_entry(metadata, header_name, coguk_columns, log_handle, column_name, count)
            elif columns is not None:
                add_header_entry(metadata, header_name, columns, log_handle, column_name, count)
            else:
                sys.exit("Must use either --gisaid or --cog_uk flag or specify columns using --columns")

            if count > 0:
                record.id += "_" + str(count)
                record.description = ""
            if record.id != '':
                SeqIO.write(record, out_fasta, "fasta-2line")

    if len(found_headers) != len(metadata[column_name].unique().tolist()):
        log_handle.write("Warning: there were %i entries in input fasta, but only %i unique headers have been added to "
                         "metadata" %(len(found_headers),len(metadata[column_name].unique().tolist())))
    metadata.to_csv(output_metadata, index=False)

    if log_handle:
        log_handle.close()
