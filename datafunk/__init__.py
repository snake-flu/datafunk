from pkg_resources import get_distribution

try:
    __version__ = get_distribution("datafunk").version
except:
    __version__ = "local"

__all__ = ["remove_dat_junk","repair_names","clean_names","remove_fasta","merge_fasta","get_ID_from_json_dict"]

from datafunk import *
