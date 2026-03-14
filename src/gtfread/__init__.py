"""Public package API for gtfread."""

try:
    from ._parser import parse_chunk_columns
except ImportError:
    from ._fallback import parse_chunk_columns

__all__ = ["parse_chunk_columns"]
