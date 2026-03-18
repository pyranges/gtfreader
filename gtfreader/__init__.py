"""Public package API for gtfreader."""

try:
    from ._parser import parse_chunk_columns
except ImportError:
    def parse_chunk_columns(*_args, **_kwargs):
        msg = "gtfreader.parse_chunk_columns requires the compiled extension. Use read_gtf_python or read_gtf_full_python for the pure Python fallback."
        raise ImportError(msg)

from .readers import (
    find_first_data_line_index,
    parse_kv_fields,
    read_gtf,
    read_gtf_full,
    read_gtf_full_python,
    read_gtf_python,
    read_gtf_restricted,
    read_gtf_restricted_python,
    to_rows,
    to_rows_keep_duplicates,
)

__all__ = [
    "find_first_data_line_index",
    "parse_chunk_columns",
    "parse_kv_fields",
    "read_gtf",
    "read_gtf_full",
    "read_gtf_full_python",
    "read_gtf_python",
    "read_gtf_restricted",
    "read_gtf_restricted_python",
    "to_rows",
    "to_rows_keep_duplicates",
]
