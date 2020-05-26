from datafunk.anonymize_microreact import *

def run(options):
        anonymize_microreact(metadata_in = options.metadata_in,
                             tree_in = options.tree_in,
                             metadata_out = options.metadata_out,
                             tree_out = options.tree_out)
