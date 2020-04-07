from pkg_resources import get_distribution

try:
    __version__ = get_distribution("datafunk").version
except:
    __version__ = "local"

__all__ = ["remove_dat_junk","repair_names","clean_names","remove_fasta","merge_fasta",
           "filter_low_coverage", "process_gisaid_sequence_data", "sam_2_fasta",
           "gisaid_json_2_metadata"]

from datafunk import *
