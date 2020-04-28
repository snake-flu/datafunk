from datafunk.mask import *

def run(options):
    mask(fasta_in = options.fasta_in,
         fasta_out = options.fasta_out,
         mask_file = options.mask_file)
