from datafunk.extract_unannotated_seqs import *

def run(options):
    extract_unannotated_seqs(fasta_in = options.fasta_in,
                       metadata_in = options.metadata_in,
                        fasta_out = options.fasta_out,
                        index_column = options.index_column,
                        null_column = options.null_column)
