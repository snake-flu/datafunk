from datafunk.exclude_uk_seqs import *

def run(options):
    exclude_UK_seqs(input = options.fasta_in,
                    output = options.fasta_out)
