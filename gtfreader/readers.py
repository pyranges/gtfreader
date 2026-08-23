"""Serial GTF readers."""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

import numpy as np
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

HASH = ord("#")


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
    return _apply_categorical_dtypes(result)


def _apply_categorical_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    """Give the five GTF categorical columns one dtype, whatever the file.

    Without this the dtype is an accident of chunking: `pd.concat` demotes a
    categorical to object when the chunks carry different categories, so on a
    coordinate-sorted GTF `Chromosome` came back object while `Feature`, whose
    values all appear in every chunk, came back category. That also made the
    result depend on `chunksize`, which is a tuning knob, not a contract.

    Category *order* matters as much as the dtype: pandas groups and sorts a
    categorical by it. `astype` sorts, and pyarrow's dictionary is in order of
    first appearance, so the two readers would otherwise disagree on the order
    of `groupby` output for the same file.
    """
    for column in GTF_DTYPES:
        if column not in df.columns:
            continue
        values = df[column]
        if isinstance(values.dtype, pd.CategoricalDtype):
            categories = values.cat.categories
            if not categories.is_monotonic_increasing:
                df[column] = values.cat.reorder_categories(categories.sort_values())
        else:
            df[column] = values.astype("category")
    return df


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
        comment="#",
        chunksize=chunksize,
        skiprows=skiprows,
        nrows=nrows,
    )


def _pyarrow_csv():
    """Return ``(pyarrow, pyarrow.csv)``, or None when pyarrow is not installed.

    The single import site for the fast path. Tests patch this to force the
    pandas reader, so there is exactly one place to disable.
    """
    try:
        import pyarrow as pa
        from pyarrow import csv as pacsv
    except ImportError:
        return None
    return pa, pacsv


def _handle_invalid_row(row):
    """Skip a whole-line comment; refuse anything else.

    pandas drops ``#`` lines because of ``comment="#"``; to pyarrow they are
    rows with the wrong number of fields. Skipping them here reproduces pandas.
    A ragged row that is *not* a comment -- a ``##FASTA`` section's sequence
    lines, say -- is a file this path does not model, so it errors and the
    caller falls back to pandas, which has its own answer for it.
    """
    return "skip" if row.text.lstrip().startswith("#") else "error"


def _contains_hash(pa, column) -> bool:
    """True when the byte ``#`` occurs anywhere in a text column.

    ``comment="#"`` truncates a line at a ``#`` in *any* position, not just at
    the start, so a file with one inside a field is read differently by the two
    parsers. Rather than model that, detect it and let pandas have the file.

    The answer only needs the raw character buffer, so this runs at memory
    speed -- about 0.03 s over the 764 MB attribute column of a 2 x 10^6-row
    GTF, against 0.7 s for ``pyarrow.compute.match_substring``. The comparison
    allocates a mask the size of one chunk, not of the column, which is why the
    chunks are read as they come rather than combined first.
    """
    for chunk in column.chunks:
        # A dictionary column keeps its text in the (tiny) dictionary.
        values = chunk.dictionary if pa.types.is_dictionary(chunk.type) else chunk
        if not (pa.types.is_string(values.type) or pa.types.is_large_string(values.type)):
            continue
        data = values.buffers()[2]
        if data is None:
            continue
        if (np.frombuffer(data, dtype=np.uint8) == HASH).any():
            return True
    return False


def _read_gtf_arrow(path: Path, *, skiprows: int, nrows: int | None):
    """Parse the nine fixed GTF columns with ``pyarrow.csv``, or return None.

    pandas' C parser is single-threaded and is about 37% of the cost of reading
    a GTF -- the attribute expansion is the rest. ``pyarrow.csv`` parses on
    every core: 1.69 s against 0.10 s for 10^6 rows of gencode on twelve cores.

    Returns None whenever pandas should read the file instead: pyarrow missing,
    a row limit (``pyarrow.csv`` has none, so reading everything to discard most
    of it would be slower), a file pyarrow cannot parse, or a ``#`` inside a
    field, where the two parsers genuinely disagree.
    """
    modules = _pyarrow_csv()
    if modules is None or nrows is not None:
        return None
    pa, pacsv = modules

    dictionary = pa.dictionary(pa.int32(), pa.string())
    try:
        table = pacsv.read_csv(
            path,
            read_options=pacsv.ReadOptions(column_names=GTF_NAMES, use_threads=True, skip_rows=skiprows),
            parse_options=pacsv.ParseOptions(delimiter="\t", invalid_row_handler=_handle_invalid_row),
            convert_options=pacsv.ConvertOptions(
                # Dictionary-encoding the five categorical columns up front
                # means the categorical arrives without a second pass. Start
                # and End are pinned rather than inferred because inference
                # turns a coordinate too large for int64 into a float, quietly
                # rounding it; pandas keeps the integer. Pinned, pyarrow
                # refuses the file instead and pandas gets it.
                column_types={
                    **{name: dictionary for name in GTF_DTYPES},
                    "Start": pa.int64(),
                    "End": pa.int64(),
                },
                # pandas applies its NA list to text columns too, so an empty
                # field reads back as NaN rather than "". Their default null
                # spellings are otherwise the same set.
                strings_can_be_null=True,
            ),
        )
    except (pa.ArrowInvalid, pa.ArrowNotImplementedError, UnicodeDecodeError, OSError):
        return None

    if any(_contains_hash(pa, table.column(name)) for name in table.column_names):
        return None

    # A file with no data lines is not worth modelling: pandas' empty frame
    # carries the dtypes of the chunk it never filled, which is not something
    # the Arrow schema reproduces. Hand the degenerate case back.
    if table.num_rows == 0:
        return None

    return _widen_empty_columns(pa, table)


def _widen_empty_columns(pa, table):
    """Give a wholly-empty column the float64 pandas would have inferred.

    An empty field is null on both sides, but a column that is *nothing but*
    empty fields has no type to infer: pandas calls it float64 full of NaN,
    pyarrow calls it null, and null converts to an object column of None. A GTF
    with no scores is the ordinary way to hit this.
    """
    fields = [
        field.with_type(pa.float64()) if pa.types.is_null(field.type) else field
        for field in table.schema
    ]
    if all(field.type == original.type for field, original in zip(fields, table.schema)):
        return table
    return table.cast(pa.schema(fields))


def _arrow_to_pandas(table) -> pd.DataFrame:
    """Convert to pandas using pandas' own spelling of a missing value.

    pyarrow puts `None` in an object column where pandas' reader puts `NaN`. A
    GTF whose Score is present on some rows and empty on others is the ordinary
    way to hit this. Under pandas 3 both are NA and the difference does not
    arise; under pandas 2 they are distinguishable, and pandas already warns
    that a future version will stop treating them as equal.

    `null_count` is Arrow metadata, so a column without nulls costs nothing.
    """
    df = table.to_pandas()
    for name in df.columns:
        if df[name].dtype == object and table.column(name).null_count:
            df[name] = df[name].fillna(np.nan)
    return df


def _frames_from_arrow_table(table, chunksize: int, parse_attributes, *, duplicate_attr: bool, ignore_bad: bool):
    """Expand attributes a chunk at a time, as the pandas path does.

    Slicing the Arrow table rather than converting it whole keeps the raw
    attribute strings of one chunk alive at a time. Every slice carries the
    file's whole dictionary, so the categorical columns survive the concat
    instead of decaying to object the way per-chunk pandas categoricals do.
    """
    frames = []
    for start in range(0, table.num_rows, chunksize):
        frame = _arrow_to_pandas(table.slice(start, chunksize))
        # pandas numbers its chunks continuously across the file; a slice
        # converted on its own would restart at zero and the concat would end
        # up with a repeated index.
        frame.index = pd.RangeIndex(start, start + len(frame))
        extra = parse_attributes(frame["Attribute"], duplicate_attr=duplicate_attr, ignore_bad=ignore_bad)
        frames.append(pd.concat([frame.drop(columns="Attribute"), extra], axis=1, sort=False))
    return frames


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
    table = _read_gtf_arrow(path, skiprows=skiprows, nrows=nrows)
    if table is not None:
        dfs = _frames_from_arrow_table(
            table,
            chunksize,
            parse_attributes,
            duplicate_attr=duplicate_attr,
            ignore_bad=ignore_bad,
        )
    else:
        dfs = []
        with _open_gtf_reader(path, chunksize=chunksize, skiprows=skiprows, nrows=nrows) as df_iter:
            for df in df_iter:
                extra = parse_attributes(df["Attribute"], duplicate_attr=duplicate_attr, ignore_bad=ignore_bad)
                dfs.append(pd.concat([df.drop(columns="Attribute"), extra], axis=1, sort=False))

    if not dfs:
        return pd.DataFrame(columns=GTF_NAMES[:-1])

    return _finalize_gtf_frame(pd.concat(dfs, sort=False))


def read_gtf(
    f: str | Path,
    /,
    *,
    nrows: int | None = None,
    duplicate_attr: bool = False,
    ignore_bad: bool = False,
) -> pd.DataFrame:
    """Read a GTF file using the compiled parser path when available."""
    path = Path(f)
    skiprows = find_first_data_line_index(path)
    return read_gtf_full(
        path,
        nrows=nrows,
        skiprows=skiprows,
        duplicate_attr=duplicate_attr,
        ignore_bad=ignore_bad,
    )


def read_gtf_python(
    f: str | Path,
    /,
    *,
    nrows: int | None = None,
    duplicate_attr: bool = False,
    ignore_bad: bool = False,
) -> pd.DataFrame:
    """Read a GTF file using the pure Python attribute parser."""
    path = Path(f)
    skiprows = find_first_data_line_index(path)
    return read_gtf_full_python(
        path,
        nrows=nrows,
        skiprows=skiprows,
        duplicate_attr=duplicate_attr,
        ignore_bad=ignore_bad,
    )


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


__all__ = [
    "find_first_data_line_index",
    "parse_kv_fields",
    "read_gtf",
    "read_gtf_full",
    "read_gtf_full_python",
    "read_gtf_python",
    "to_rows",
    "to_rows_keep_duplicates",
]
