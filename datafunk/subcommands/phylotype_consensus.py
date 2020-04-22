from datafunk.phylotype_consensus import *

def run(options):
    align_by_phylotype(options.input_fasta,options.clade_file,options.input_metadata,options.output_folder)