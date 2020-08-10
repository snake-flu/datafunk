from datafunk.bootstrap import *

def run(options):
    bootstrap(fasta_in = options.fasta_in,
              n = options.n,
              output_prefix = options.output_prefix)
