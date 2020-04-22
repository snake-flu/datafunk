from Bio import SeqIO
import argparse

def pad_alignment(alignment, leftpad, rightpad, output):
    if output:
        out = open(output, 'w')
    else:
        out = sys.stdout

    with open(alignment, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            id = record.id
            seq = 'N' * leftpad + str(record.seq) + 'N' * rightpad
            out.write('>' + id + '\n')
            out.write(seq + '\n')

    if output:
        out.close()

    pass
