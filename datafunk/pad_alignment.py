from Bio import SeqIO
import argparse, sys

def pad_alignment(alignment, leftpad, rightpad, output):
    if output:
        out = open(output, 'w')
    else:
        out = sys.stdout

    with open(alignment, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            id = record.id
            seq = 'N' * int(leftpad) + str(record.seq) + 'N' * int(rightpad)
            out.write('>' + id + '\n')
            out.write(seq + '\n')

    if output:
        out.close()

    pass
