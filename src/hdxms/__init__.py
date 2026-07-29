"""Reusable HDX-MS analysis and structure-coloring utilities."""
from .consensus import consensus_from_csv, consensus_from_dataframe
from .deduplicate import remove_duplicates

__all__ = ["consensus_from_csv", "consensus_from_dataframe", "remove_duplicates"]
__version__ = "1.2.0"
