from datafunk.repair_names import *

def run(options):
    fix_names(options.fasta, options.tree, options.out)