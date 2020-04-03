from datafunk.split_by_phylotype import *

def run(options):
    phylo_dictionary = import_file(options.input_metafile,options.input_fasta)
    split_phylotype(options.output_folder,phylo_dictionary,options.clade_threshold)