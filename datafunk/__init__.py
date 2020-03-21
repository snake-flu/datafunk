from pkg_resources import get_distribution

try:
    __version__ = get_distribution("datafunk").version
except:
    __version__ = "local"

__all__ = ["remove_dat_junk"]

from datafunk import *
