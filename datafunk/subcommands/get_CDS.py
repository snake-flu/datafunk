from datafunk.get_CDS import *

def run(options):
    get_CDS(fasta_in = options.fasta_in,
            fasta_out = options.fasta_out,
            translate = options.translate)
