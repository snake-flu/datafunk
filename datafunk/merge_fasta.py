from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import glob
import csv

def merge_fasta(input_folder, metafile, output_file="merged_file.fasta"):
    metadata_dictionary = {}
    sequence_dictionary = {}
    merged_file = open(output_file,"w")
    log_file = open(output_file+".log","w")
    with open(metafile) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            metadata_dictionary[row[0]] = row[1:]

    for fasta_file in glob.glob(input_folder+"/*.fasta"):
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in metadata_dictionary.keys() and record.id not in sequence_dictionary.keys():
                sequence_dictionary[str(record.id).rstrip()] = record.seq
            elif record.id not in metadata_dictionary.keys():
                log_file.write(record.id + " sequence is not in metadata file or the name is wrong (file " + fasta_file + ")\n")
            elif record.id in sequence_dictionary.keys():
                log_file.write(record.id + " is a duplicate (file " + fasta_file + ")\n")

    for key, value in sequence_dictionary.items():
        records = SeqRecord(value, key, description= '')
        SeqIO.write(records, merged_file, "fasta")

    merged_file.close()
    log_file.close()