from datafunk.pad_alignment import *

def run(options):
    if options.stdout:
        fasta_out = None
    else:
        fasta_out = options.fasta_out

    pad_alignment(alignment = options.fasta_in,
                  leftpad = options.left_pad,
                  rightpad = options.right_pad,
                  output = fasta_out)
