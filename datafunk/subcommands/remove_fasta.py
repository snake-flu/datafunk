from datafunk.remove_fasta import *

def run(options):
    filter_dictionary = filter_list(options.filter_file)
    remove_fasta(options.input_fasta,filter_dictionary,options.output_file)
