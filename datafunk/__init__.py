from pkg_resources import get_distribution

try:
    __version__ = get_distribution("datafunk").version
except:
    __version__ = "local"

__all__ = ["repair_names","clean_names","remove_fasta","merge_fasta",
           "filter_fasta_by_covg_and_length", "process_gisaid_sequence_data", "sam_2_fasta",
           "phylotype_consensus", "gisaid_json_2_metadata", "set_uniform_header", "add_epi_week",
           "process_gisaid_data", "pad_alignment", "exclude_uk_seqs", "get_CDS",
           "distance_to_root", "mask", "curate_lineages", "snp_finder", "add_header_column",
           "anonymize_microreact", "extract_unannotated_seqs", "del_finder", "AA_finder"]

from datafunk import *
