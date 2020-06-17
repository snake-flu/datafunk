from datafunk.del_finder import *

def run(options):
    del_finder(fasta_in = options.fasta_in,
               fasta_out = options.fasta_out,
               del_file = options.del_file,
               genotypes_file = options.genotypes_file,
               append_snp = options.append_snp)
