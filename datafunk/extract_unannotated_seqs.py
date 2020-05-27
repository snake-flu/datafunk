import pandas as pd
from Bio import SeqIO

def extract_unannotated_seqs(fasta_in, fasta_out, metadata_in, index_column, null_column):
    fasta_in = SeqIO.index(str(fasta_in), "fasta")
    df = pd.read_csv(metadata_in)

    with open(str(fasta_out), 'w') as fasta_out:
        for i,row in df.iterrows():
            if pd.isnull(row[null_column]):
                sequence_name = row[index_column]
                if sequence_name in fasta_in:
                    record = fasta_in[sequence_name]
                    fasta_out.write('>' + record.id + '\n')
                    fasta_out.write(str(record.seq) + '\n')
