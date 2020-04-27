from datafunk.distance_to_root import *

def run(options):
    distance_to_root(fasta_file = options.fasta_in,
                     metadata_file = options.metadata_in,
                     plot = options.plot)
