from datafunk.snp_finder import *

def run(options):
    read_alignment_and_get_snps(options.a, options.snp, options.o)