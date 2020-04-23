from Bio import SeqIO
import sys

def exclude_UK_seqs(input, output):
    out_handle = open(output, 'w')
    with open(input, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            id = record.id
            seq = str(record.seq)
            if id.split('/')[0] in ['England', 'Wales', 'Scotland', 'Northern_Ireland']:
                continue
            else:
                out_handle.write('>' + id + '\n')
                out_handle.write(seq + '\n')

    out_handle.close()
