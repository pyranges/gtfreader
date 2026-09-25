"""The pyarrow path must return exactly what the pandas path returns.

pandas is the reference implementation: whatever it returns today is the
behaviour users have. Every test here runs twice, once with the fast path
disabled, and the parity tests compare the two frames strictly -- dtypes and
category *order* included, because pandas groups and sorts a categorical by its
category order, so getting that wrong silently reorders `groupby` output.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from gtfreader import read_gtf, read_gtf_full
from gtfreader import readers

DATA_LINE = 'chr1\thavana\tgene\t11869\t14409\t.\t+\t.\tgene_id "G1"; gene_name "DDX11L1";\n'

CORPUS = {
    "bare": DATA_LINE,
    "leading_comments": "#!genome-build GRCh38\n##provider: GENCODE\n" + DATA_LINE,
    "comments_in_body": (
        "# header\n"
        + DATA_LINE
        + "### separator\n"
        + 'chr2\thavana\texon\t20\t30\t.\t-\t.\tgene_id "G2"; exon_number "1";\n'
        + "##sequence-region chr3 1 100\n"
        + 'chr3\tensembl\tCDS\t40\t50\t.\t+\t0\tgene_id "G3";\n'
    ),
    "unsorted_chromosomes": (
        'chr9\thavana\tgene\t1\t10\t.\t+\t.\tgene_id "G9";\n'
        'chr1\thavana\tgene\t2\t20\t.\t-\t.\tgene_id "G1";\n'
        'chr22\thavana\tgene\t3\t30\t.\t+\t.\tgene_id "G22";\n'
        'chr2\thavana\tgene\t4\t40\t.\t-\t.\tgene_id "G2";\n'
    ),
    "many_rows": "".join(
        f'chr{(i % 7) + 1}\thavana\tgene\t{i + 1}\t{i + 100}\t.\t+\t.\tgene_id "G{i}"; level "{i % 3}";\n'
        for i in range(500)
    ),
    "ragged_attributes": (
        'chr1\thavana\tgene\t1\t10\t.\t+\t.\tgene_id "G1";\n'
        'chr1\thavana\ttranscript\t1\t10\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; tag "basic";\n'
        'chr1\thavana\texon\t1\t10\t.\t+\t.\texon_id "E1";\n'
    ),
    "empty_attribute": 'chr1\thavana\tgene\t11869\t14409\t.\t+\t.\t\n',
    "numeric_score": 'chr1\thavana\tgene\t1\t10\t42\t+\t.\tgene_id "G1";\n',
    "no_trailing_newline": DATA_LINE.rstrip("\n"),
    "blank_lines": DATA_LINE + "\n\n" + 'chr2\thavana\tgene\t5\t9\t.\t-\t.\tgene_id "G2";\n',
    "duplicate_attributes": (
        'chr1\thavana\texon\t1\t10\t.\t+\t.\tgene_id "G1"; tag "CCDS"; tag "basic";\n'
    ),
    # Every case below was found by comparing the two readers on files the
    # handwritten cases above did not reach. Each one diverged at least once.
    "empty_score_column": 'chr1\thavana\tgene\t1\t10\t\t+\t.\tgene_id "G1";\n' * 3,
    "empty_attribute_column": 'chr1\thavana\tgene\t1\t10\t.\t+\t.\t\n' * 3,
    "only_comments": "# a\n## b\n",
    "empty_file": "",
    "score_all_numeric": 'chr1\thavana\tgene\t1\t10\t42\t+\t.\tgene_id "G1";\n' * 3,
    "score_partly_empty": (
        'chr1\thavana\tgene\t1\t10\t.\t+\t.\tgene_id "G1";\n'
        'chr1\thavana\tgene\t2\t20\t\t+\t.\tgene_id "G2";\n'
    ),
    "score_mixed": DATA_LINE + 'chr1\thavana\tgene\t1\t10\t42\t+\t.\tgene_id "G2";\n',
    "coordinate_too_large_for_int64": 'chr1\thavana\tgene\t1\t99999999999999999999\t.\t+\t.\tgene_id "G1";\n',
    "na_spelling_in_score": 'chr1\thavana\tgene\t1\t10\tNA\t+\t.\tgene_id "G1";\n' * 3,
    "chromosome_named_na": 'NA\thavana\tgene\t1\t10\t.\t+\t.\tgene_id "G1";\n' * 3,
    "utf8_attribute": 'chr1\thavana\tgene\t1\t10\t.\t+\t.\tgene_id "éèü";\n',
    "attribute_starting_with_a_quote": 'chr1\thavana\tgene\t1\t10\t.\t+\t.\t"quoted"; gene_id "G1";\n',
    "unbalanced_quote": 'chr1\thavana\tgene\t1\t10\t.\t+\t.\tgene_id "a"b"; x "y";\n',
    "too_many_fields": 'chr1\thavana\tgene\t1\t10\t.\t+\t.\tgene_id "G1";\textra\n',
    "too_few_fields": "chr1\thavana\tgene\t1\t10\t.\t+\t.\n",
    "crlf": 'chr1\thavana\tgene\t1\t10\t.\t+\t.\tgene_id "G1";\r\n' * 3,
}


@pytest.fixture(params=["pandas", "pyarrow"])
def mode(request, monkeypatch) -> str:
    """Run the test body once per reader.

    `_pyarrow_csv` is the single import site for the fast path, so patching it
    is the whole of "pretend pyarrow is not installed".
    """
    if request.param == "pandas":
        monkeypatch.setattr(readers, "_pyarrow_csv", lambda: None)
    elif readers._pyarrow_csv() is None:
        pytest.skip("pyarrow is not installed")
    return request.param


@pytest.fixture
def requires_pyarrow():
    if readers._pyarrow_csv() is None:
        pytest.skip("pyarrow is not installed")


def write_gtf(tmp_path: Path, contents: str, *, name: str = "test.gtf", compress: bool = False) -> Path:
    path = tmp_path / (name + ".gz" if compress else name)
    if compress:
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(contents)
    else:
        path.write_text(contents, encoding="utf-8")
    return path


def read_with_pandas(path: Path, monkeypatch, **kwargs) -> pd.DataFrame:
    monkeypatch.setattr(readers, "_pyarrow_csv", lambda: None)
    try:
        return read_gtf(path, **kwargs)
    finally:
        monkeypatch.undo()


def read_both(path: Path, monkeypatch, **kwargs) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(pandas frame, pyarrow frame) for the same file."""
    return read_with_pandas(path, monkeypatch, **kwargs), read_gtf(path, **kwargs)


# --------------------------------------------------------------------------
# Parity
# --------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(CORPUS))
@pytest.mark.parametrize("compress", [False, True], ids=["plain", "gzip"])
def test_both_readers_agree(tmp_path, monkeypatch, requires_pyarrow, name, compress):
    path = write_gtf(tmp_path, CORPUS[name], compress=compress)
    from_pandas, from_arrow = read_both(path, monkeypatch)
    assert_frame_equal(from_pandas, from_arrow)


@pytest.mark.parametrize("duplicate_attr", [False, True], ids=["last", "duplicates"])
@pytest.mark.parametrize("ignore_bad", [False, True], ids=["strict", "lenient"])
def test_both_readers_agree_on_attribute_options(
    tmp_path, monkeypatch, requires_pyarrow, duplicate_attr, ignore_bad
):
    path = write_gtf(tmp_path, CORPUS["duplicate_attributes"] + CORPUS["ragged_attributes"])
    from_pandas, from_arrow = read_both(
        path, monkeypatch, duplicate_attr=duplicate_attr, ignore_bad=ignore_bad
    )
    assert_frame_equal(from_pandas, from_arrow)


@pytest.mark.parametrize("chunksize", [1, 2, 7, 100_000])
def test_result_does_not_depend_on_chunksize(tmp_path, mode, chunksize):
    """`chunksize` is a tuning knob; it must not reach the returned frame.

    It used to: `pd.concat` demotes a categorical to object when the chunks
    carry different categories, so the dtype of Chromosome depended on how the
    rows happened to fall across chunk boundaries.
    """
    path = write_gtf(tmp_path, CORPUS["unsorted_chromosomes"])
    chunked = read_gtf_full(path, chunksize=chunksize)
    whole = read_gtf_full(path, chunksize=1_000_000)
    assert_frame_equal(chunked, whole)


# --------------------------------------------------------------------------
# The dtype contract
# --------------------------------------------------------------------------


@pytest.mark.parametrize("chunksize", [1, 2, 100_000])
def test_categorical_columns_are_categorical(tmp_path, mode, chunksize):
    path = write_gtf(tmp_path, CORPUS["unsorted_chromosomes"])
    frame = read_gtf_full(path, chunksize=chunksize)
    for column in ("Chromosome", "Source", "Feature", "Strand", "Frame"):
        assert isinstance(frame[column].dtype, pd.CategoricalDtype), f"{column} at chunksize={chunksize}"


def test_categories_are_sorted(tmp_path, mode):
    """pyarrow's dictionary is in order of first appearance; pandas sorts.

    Category order is not cosmetic -- `groupby` and `sort_values` on a
    categorical follow it -- so the two readers have to agree on it.
    """
    path = write_gtf(tmp_path, CORPUS["unsorted_chromosomes"])
    categories = list(read_gtf(path)["Chromosome"].cat.categories)
    assert categories == sorted(categories)
    assert categories == ["chr1", "chr2", "chr22", "chr9"]


def test_index_is_continuous_across_chunks(tmp_path, mode):
    path = write_gtf(tmp_path, CORPUS["many_rows"])
    frame = read_gtf_full(path, chunksize=7)
    assert list(frame.index) == list(range(len(frame)))


# --------------------------------------------------------------------------
# When the fast path must decline
# --------------------------------------------------------------------------


def test_fast_path_is_really_taken(tmp_path, requires_pyarrow):
    """Without this the pyarrow half of every test above could be a no-op."""
    path = write_gtf(tmp_path, CORPUS["comments_in_body"])
    assert readers._read_gtf_arrow(path, skiprows=0, nrows=None) is not None


@pytest.mark.parametrize(
    ("name", "contents", "kwargs"),
    [
        (
            "hash inside an attribute",
            'chr1\thavana\tgene\t1\t10\t.\t+\t.\tgene_id "G1"; note "grade #2";\n',
            {},
        ),
        (
            "hash inside a comment line with nine fields",
            "#chr1\tx\ty\t1\t2\t.\t+\t.\tz\n" + DATA_LINE,
            {},
        ),
        (
            "fasta section",
            DATA_LINE + "##FASTA\n>chr1\nACGTACGT\n",
            {},
        ),
        ("row limit", DATA_LINE * 10, {"nrows": 3}),
    ],
)
def test_fast_path_declines(tmp_path, requires_pyarrow, name, contents, kwargs):
    path = write_gtf(tmp_path, contents)
    assert readers._read_gtf_arrow(path, skiprows=0, nrows=kwargs.get("nrows")) is None, name


@pytest.mark.parametrize(
    ("name", "contents"),
    [
        ("hash inside an attribute", 'chr1\thavana\tgene\t1\t10\t.\t+\t.\tgene_id "G1"; note "grade #2";\n'),
        ("fasta section", DATA_LINE + "##FASTA\n>chr1\nACGTACGT\n"),
    ],
)
def test_declined_files_still_read(tmp_path, monkeypatch, requires_pyarrow, name, contents):
    """Declining means pandas reads it, not that the read fails."""
    path = write_gtf(tmp_path, contents)
    from_pandas, from_arrow = read_both(path, monkeypatch)
    assert_frame_equal(from_pandas, from_arrow)


def test_readers_fail_the_same_way(tmp_path, monkeypatch, requires_pyarrow):
    """A file neither reader can make sense of must fail identically."""
    path = write_gtf(tmp_path, 'chr1\thavana\tgene\t.\t10\t.\t+\t.\tgene_id "G1";\n')
    with pytest.raises(TypeError) as from_pandas:
        read_with_pandas(path, monkeypatch)
    with pytest.raises(TypeError) as from_arrow:
        read_gtf(path)
    assert type(from_pandas.value) is type(from_arrow.value)


def test_nulls_use_pandas_spelling(tmp_path, monkeypatch, requires_pyarrow):
    """pyarrow puts None in an object column where pandas' reader puts NaN.

    Indistinguishable under pandas 3, where both are NA, and distinguishable
    under pandas 2 -- which also warns that a future version will stop treating
    them as equal.
    """
    path = write_gtf(tmp_path, CORPUS["score_partly_empty"])
    from_pandas, from_arrow = read_both(path, monkeypatch)
    assert from_pandas["Score"].isna().sum() > 0
    assert {type(v) for v in from_pandas["Score"][from_pandas["Score"].isna()]} == {
        type(v) for v in from_arrow["Score"][from_arrow["Score"].isna()]
    }


def test_wholly_empty_column_is_not_object(tmp_path, mode):
    """pandas infers float64 for a column of nothing but empty fields.

    pyarrow infers its `null` type, which converts to an object column of None
    -- same information, different dtype, and dtype is the contract.
    """
    path = write_gtf(tmp_path, CORPUS["empty_score_column"])
    assert read_gtf(path)["Score"].dtype == "float64"


def test_nrows_agrees_between_readers(tmp_path, monkeypatch, requires_pyarrow):
    path = write_gtf(tmp_path, CORPUS["many_rows"])
    from_pandas, from_arrow = read_both(path, monkeypatch, nrows=13)
    assert len(from_arrow) == 13
    assert_frame_equal(from_pandas, from_arrow)


def test_reader_works_without_pyarrow(tmp_path, monkeypatch):
    """The import failing is the supported state, not an error path."""
    import builtins

    real_import = builtins.__import__

    def no_pyarrow(name, *args, **kwargs):
        if name.startswith("pyarrow"):
            msg = "No module named 'pyarrow'"
            raise ImportError(msg)
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_pyarrow)
    assert readers._pyarrow_csv() is None

    path = write_gtf(tmp_path, CORPUS["comments_in_body"])
    frame = read_gtf(path)
    assert list(frame["gene_id"]) == ["G1", "G2", "G3"]


def test_comment_lines_are_dropped_not_parsed(tmp_path, mode):
    path = write_gtf(tmp_path, CORPUS["comments_in_body"])
    frame = read_gtf(path)
    assert list(frame["Chromosome"]) == ["chr1", "chr2", "chr3"]
    assert not any(str(value).startswith("#") for value in frame["Chromosome"])


def test_start_is_converted_to_zero_based(tmp_path, mode):
    path = write_gtf(tmp_path, DATA_LINE)
    assert read_gtf(path).iloc[0]["Start"] == 11868


# --------------------------------------------------------------------------
# The public switch
# --------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(CORPUS))
def test_use_pyarrow_false_matches_the_default(tmp_path, requires_pyarrow, name):
    """Forcing the pandas parser changes speed, never the frame.

    The parity tests above reach the slow path by patching a private function.
    This one goes through the documented keyword, so the switch callers
    actually have is the one under test.
    """
    path = write_gtf(tmp_path, CORPUS[name])
    assert_frame_equal(read_gtf(path), read_gtf(path, use_pyarrow=False))
    assert_frame_equal(read_gtf_full(path), read_gtf_full(path, use_pyarrow=False))


def test_use_pyarrow_false_really_declines(tmp_path, requires_pyarrow):
    path = write_gtf(tmp_path, CORPUS["bare"])
    assert readers._read_gtf_arrow(path, skiprows=0, nrows=None, use_pyarrow=False) is None
    assert readers._read_gtf_arrow(path, skiprows=0, nrows=None) is not None


def test_use_pyarrow_true_requires_pyarrow(tmp_path, monkeypatch):
    """Asking for the fast parse and silently getting the slow one is worse than an error.

    Nothing is broken for anyone who leaves the default alone: without pyarrow,
    `None` still falls back. Only an explicit `True` insists.
    """
    monkeypatch.setattr(readers, "_pyarrow_csv", lambda: None)
    path = write_gtf(tmp_path, CORPUS["bare"])
    with pytest.raises(ImportError, match="pyarrow is not installed"):
        read_gtf(path, use_pyarrow=True)
    with pytest.raises(ImportError, match="pyarrow is not installed"):
        read_gtf_full(path, use_pyarrow=True)
    # The default is unaffected -- this is the install that must keep working.
    assert_frame_equal(read_gtf(path), read_gtf(path, use_pyarrow=False))


def test_use_pyarrow_true_with_nrows_still_reads(tmp_path, requires_pyarrow):
    """`nrows` declines the fast path for speed, not for lack of pyarrow.

    Raising there would make `use_pyarrow=True` mean two different things
    depending on an unrelated argument.
    """
    path = write_gtf(tmp_path, CORPUS["bare"])
    assert_frame_equal(read_gtf(path, use_pyarrow=True, nrows=1), read_gtf(path, nrows=1))
