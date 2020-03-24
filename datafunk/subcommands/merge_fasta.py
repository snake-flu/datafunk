from datafunk.merge_fasta import *

def run(options):
    merge_fasta(options.folder, options.input_metafile, options.output_file)
