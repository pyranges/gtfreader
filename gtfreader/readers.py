"""Serial GTF readers."""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

import pandas as pd

try:
    from ._parser import parse_chunk_columns as _parse_chunk_columns_compiled
except ImportError:
    _parse_chunk_columns_compiled = None

LOGGER = logging.getLogger(__name__)

GTF_DTYPES = {
    "Chromosome": "category",
    "Source": "category",
    "Feature": "category",
    "Strand": "category",
    "Frame": "category",
}
GTF_NAMES = ["Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attribute"]
RESTRICTED_ATTRIBUTE_COLUMNS = ["gene_id", "transcript_id", "exon_number", "exon_id"]


def find_first_data_line_index(file_path: str | Path) -> int:
    """Find the first non-empty line that is not a comment."""
    path = Path(file_path)
    opener = gzip.open if path.suffix == ".gz" else open

    with opener(path, "rt", encoding="utf-8") as handle:
        for index, line in enumerate(handle):
            stripped = line.lstrip()
            if stripped and not stripped.startswith("#"):
                return index

    return 0


def parse_kv_fields(line: str) -> list[tuple[str, str]]:
    """Parse a GTF attribute string with quoted or unquoted values."""
    fields: list[tuple[str, str]] = []
    n = len(line)
    pos = 0

    while pos < n:
        while pos < n and line[pos] in " \t;":
            pos += 1
        if pos >= n:
            break

        key_start = pos
        while pos < n and line[pos] not in " \t;":
            pos += 1
        key = line[key_start:pos]

        while pos < n and line[pos] in " \t":
            pos += 1
        if pos >= n:
            break

        if line[pos] == '"':
            pos += 1
            value_start = pos
            while pos < n and line[pos] != '"':
                pos += 1
            value = line[value_start:pos]
            if pos < n and line[pos] == '"':
                pos += 1
        else:
            value_start = pos
            while pos < n and line[pos] != ';':
                pos += 1
            value = line[value_start:pos].strip()

        fields.append((key, value))

        while pos < n and line[pos] != ';':
            pos += 1
        if pos < n and line[pos] == ';':
            pos += 1

    return fields


def to_rows(attribute_column: pd.Series, *, ignore_bad: bool = False) -> pd.DataFrame:
    """Parse a GTF attribute column into a dataframe of attribute columns."""
    attribute_column = _normalize_attribute_series(attribute_column)
    rowdicts: list[dict[str, str]] = []
    line = ""
    try:
        for line in attribute_column:
            rowdicts.append(dict(parse_kv_fields(line)))
    except ValueError:
        if not ignore_bad:
            LOGGER.exception(
                "The following line is not parseable as GTF:\n%s\n\nTo ignore bad lines use ignore_bad=True.",
                line,
            )
            raise

    return pd.DataFrame.from_records(rowdicts, index=attribute_column.index)


def to_rows_keep_duplicates(attribute_column: pd.Series, *, ignore_bad: bool = False) -> pd.DataFrame:
    """Parse a GTF attribute column and keep duplicate attributes as comma-joined values."""
    attribute_column = _normalize_attribute_series(attribute_column)
    rowdicts: list[dict[str, str]] = []
    line = ""
    try:
        for line in attribute_column:
            rowdict: dict[str, list[str]] = {}
            for key, value in parse_kv_fields(line):
                rowdict.setdefault(key, []).append(value)
            rowdicts.append({key: ",".join(values) for key, values in rowdict.items()})
    except ValueError:
        if not ignore_bad:
            LOGGER.exception(
                "The following line is not parseable as GTF:\n%s\n\nTo ignore bad lines use ignore_bad=True.",
                line,
            )
            raise

    return pd.DataFrame.from_records(rowdicts, index=attribute_column.index)


def _normalize_attribute_series(attribute_column: pd.Series) -> pd.Series:
    return attribute_column.where(attribute_column.notna(), "").astype(str)


def _finalize_gtf_frame(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df

    result = df.copy()
    result.loc[:, "Start"] = result["Start"] - 1
    return result


def _resolve_chunksize(chunksize: int, chunk_size: int | None) -> int:
    if chunk_size is None:
        return chunksize
    if chunk_size <= 0:
        msg = "chunk_size must be greater than 0."
        raise ValueError(msg)
    return chunk_size


def _open_gtf_reader(
    path: Path,
    *,
    chunksize: int,
    skiprows: int,
    nrows: int | None,
):
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=GTF_NAMES,
        dtype=GTF_DTYPES,
        chunksize=chunksize,
        skiprows=skiprows,
        nrows=nrows,
    )


def _parse_attributes_compiled(
    attribute_column: pd.Series,
    *,
    duplicate_attr: bool = False,
    ignore_bad: bool = False,
) -> pd.DataFrame:
    if _parse_chunk_columns_compiled is None:
        return _parse_attributes_python(
            attribute_column,
            duplicate_attr=duplicate_attr,
            ignore_bad=ignore_bad,
        )

    attribute_column = _normalize_attribute_series(attribute_column)
    if duplicate_attr:
        return to_rows_keep_duplicates(attribute_column, ignore_bad=ignore_bad)
    if ignore_bad:
        return to_rows(attribute_column, ignore_bad=ignore_bad)
    return pd.DataFrame(
        _parse_chunk_columns_compiled(attribute_column.to_numpy(copy=False)),
        index=attribute_column.index,
    )


def _parse_attributes_python(
    attribute_column: pd.Series,
    *,
    duplicate_attr: bool = False,
    ignore_bad: bool = False,
) -> pd.DataFrame:
    attribute_column = _normalize_attribute_series(attribute_column)
    if duplicate_attr:
        return to_rows_keep_duplicates(attribute_column, ignore_bad=ignore_bad)
    return to_rows(attribute_column, ignore_bad=ignore_bad)


def _parse_restricted_attributes_compiled(attribute_column: pd.Series) -> pd.DataFrame:
    if _parse_chunk_columns_compiled is None:
        return _parse_restricted_attributes_python(attribute_column)

    attribute_column = _normalize_attribute_series(attribute_column)
    columns = _parse_chunk_columns_compiled(attribute_column.to_numpy(copy=False))
    return pd.DataFrame(
        {name: columns[name] for name in RESTRICTED_ATTRIBUTE_COLUMNS},
        index=attribute_column.index,
    )


def _parse_restricted_attributes_python(attribute_column: pd.Series) -> pd.DataFrame:
    return to_rows(_normalize_attribute_series(attribute_column)).reindex(columns=RESTRICTED_ATTRIBUTE_COLUMNS)


def _read_gtf_full(
    path: Path,
    *,
    nrows: int | None,
    skiprows: int,
    chunksize: int,
    duplicate_attr: bool,
    ignore_bad: bool,
    parse_attributes,
) -> pd.DataFrame:
    dfs: list[pd.DataFrame] = []
    with _open_gtf_reader(path, chunksize=chunksize, skiprows=skiprows, nrows=nrows) as df_iter:
        for df in df_iter:
            extra = parse_attributes(df["Attribute"], duplicate_attr=duplicate_attr, ignore_bad=ignore_bad)
            dfs.append(pd.concat([df.drop(columns="Attribute"), extra], axis=1, sort=False))

    if not dfs:
        return pd.DataFrame(columns=GTF_NAMES[:-1])

    return _finalize_gtf_frame(pd.concat(dfs, sort=False))


def _read_gtf_restricted(
    path: Path,
    *,
    skiprows: int,
    nrows: int | None,
    chunksize: int,
    parse_attributes,
) -> pd.DataFrame:
    dfs: list[pd.DataFrame] = []
    with _open_gtf_reader(path, chunksize=chunksize, skiprows=skiprows, nrows=nrows) as df_iter:
        for df in df_iter:
            subset = parse_attributes(df["Attribute"])
            dfs.append(
                pd.concat(
                    [df[["Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame"]], subset],
                    axis=1,
                    sort=False,
                )
            )

    if not dfs:
        return pd.DataFrame(columns=GTF_NAMES[:-1] + RESTRICTED_ATTRIBUTE_COLUMNS)

    return _finalize_gtf_frame(pd.concat(dfs, sort=False))


def read_gtf(
    f: str | Path,
    /,
    *,
    nrows: int | None = None,
    full: bool = True,
    duplicate_attr: bool = False,
    ignore_bad: bool = False,
) -> pd.DataFrame:
    """Read a GTF file using the compiled parser path when available."""
    path = Path(f)
    skiprows = find_first_data_line_index(path)

    if full:
        return read_gtf_full(
            path,
            nrows=nrows,
            skiprows=skiprows,
            duplicate_attr=duplicate_attr,
            ignore_bad=ignore_bad,
        )

    return read_gtf_restricted(path, skiprows=skiprows, nrows=nrows)


def read_gtf_python(
    f: str | Path,
    /,
    *,
    nrows: int | None = None,
    full: bool = True,
    duplicate_attr: bool = False,
    ignore_bad: bool = False,
) -> pd.DataFrame:
    """Read a GTF file using the pure Python attribute parser."""
    path = Path(f)
    skiprows = find_first_data_line_index(path)

    if full:
        return read_gtf_full_python(
            path,
            nrows=nrows,
            skiprows=skiprows,
            duplicate_attr=duplicate_attr,
            ignore_bad=ignore_bad,
        )

    return read_gtf_restricted_python(path, skiprows=skiprows, nrows=nrows)


def read_gtf_full(
    f: str | Path,
    /,
    nrows: int | None = None,
    skiprows: int = 0,
    chunksize: int = int(1e5),
    *,
    chunk_size: int | None = None,
    duplicate_attr: bool = False,
    ignore_bad: bool = False,
) -> pd.DataFrame:
    """Read a GTF file and expand the attribute column using the compiled parser path."""
    path = Path(f)
    chunksize = _resolve_chunksize(chunksize, chunk_size)
    return _read_gtf_full(
        path,
        nrows=nrows,
        skiprows=skiprows,
        chunksize=chunksize,
        duplicate_attr=duplicate_attr,
        ignore_bad=ignore_bad,
        parse_attributes=_parse_attributes_compiled,
    )


def read_gtf_full_python(
    f: str | Path,
    /,
    nrows: int | None = None,
    skiprows: int = 0,
    chunksize: int = int(1e5),
    *,
    chunk_size: int | None = None,
    duplicate_attr: bool = False,
    ignore_bad: bool = False,
) -> pd.DataFrame:
    """Read a GTF file and expand the attribute column using the pure Python parser."""
    path = Path(f)
    chunksize = _resolve_chunksize(chunksize, chunk_size)
    return _read_gtf_full(
        path,
        nrows=nrows,
        skiprows=skiprows,
        chunksize=chunksize,
        duplicate_attr=duplicate_attr,
        ignore_bad=ignore_bad,
        parse_attributes=_parse_attributes_python,
    )


def read_gtf_restricted(
    f: str | Path,
    /,
    skiprows: int = 0,
    nrows: int | None = None,
    chunksize: int = int(1e5),
    *,
    chunk_size: int | None = None,
) -> pd.DataFrame:
    """Read core GTF columns plus a small compiled-parser attribute subset."""
    path = Path(f)
    chunksize = _resolve_chunksize(chunksize, chunk_size)
    return _read_gtf_restricted(
        path,
        skiprows=skiprows,
        nrows=nrows,
        chunksize=chunksize,
        parse_attributes=_parse_restricted_attributes_compiled,
    )


def read_gtf_restricted_python(
    f: str | Path,
    /,
    skiprows: int = 0,
    nrows: int | None = None,
    chunksize: int = int(1e5),
    *,
    chunk_size: int | None = None,
) -> pd.DataFrame:
    """Read core GTF columns plus a small pure-Python attribute subset."""
    path = Path(f)
    chunksize = _resolve_chunksize(chunksize, chunk_size)
    return _read_gtf_restricted(
        path,
        skiprows=skiprows,
        nrows=nrows,
        chunksize=chunksize,
        parse_attributes=_parse_restricted_attributes_python,
    )


__all__ = [
    "find_first_data_line_index",
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
