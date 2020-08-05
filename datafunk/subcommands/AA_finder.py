from datafunk.AA_finder import *

def run(options):
    AA_finder(fasta_in = options.fasta_in,
               AA_file = options.AA_file,
               genotypes_file = options.genotypes_file)
