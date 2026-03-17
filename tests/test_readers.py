from pathlib import Path

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from gtfread import (
    find_first_data_line_index,
    parse_kv_fields,
    read_gtf,
    read_gtf_full,
    read_gtf_full_python,
    read_gtf_python,
    read_gtf_restricted,
    read_gtf_restricted_python,
)


def _normalize_frame_for_compare(df):
    df = df.copy()
    for column in ["Chromosome", "Source", "Feature", "Strand", "Frame"]:
        if column in df.columns:
            df[column] = df[column].astype(str)
    return df.astype(object).where(pd.notna(df), None)


def _write_temp_gtf(tmp_path: Path, contents: str) -> Path:
    path = tmp_path / "test.gtf"
    path.write_text(contents)
    return path


@pytest.mark.parametrize(
    ("line", "expected"),
    [
        (
            'gene_id "G1"; note "contains; semicolon"; exon_number "1";',
            [("gene_id", "G1"), ("note", "contains; semicolon"), ("exon_number", "1")],
        ),
        (
            'gene_id "G2"; transcript_id "T2"; exon_number 2; level 3;',
            [("gene_id", "G2"), ("transcript_id", "T2"), ("exon_number", "2"), ("level", "3")],
        ),
    ],
)
def test_parse_kv_fields_supports_gtf_format_variants(line: str, expected: list[tuple[str, str]]):
    assert parse_kv_fields(line) == expected


def test_find_first_data_line_index_skips_comments(tmp_path: Path):
    path = _write_temp_gtf(tmp_path, "# comment\n## comment\nchr1\tsrc\tgene\t1\t2\t.\t+\t.\tgene_id \"G1\";\n")
    assert find_first_data_line_index(path) == 2


def test_read_gtf_parses_semicolons_inside_quoted_attributes(tmp_path: Path):
    path = _write_temp_gtf(
        tmp_path,
        "chr1\tRefSeq\tCDS\t20486313\t20486315\t.\t+\t0\t"
        'gene_id "SELENOO"; transcript_id "NM_001115017.5"; db_xref "GeneID:417745"; '
        'gbkey "CDS"; gene "SELENOO"; note "UGA stop codon recoded as selenocysteine; '
        'The RefSeq protein has 1 substitution compared to this genomic sequence"; '
        'product "protein adenylyltransferase SelO, mitochondrial"; '
        'protein_id "NP_001108489.5"; transl_except "(pos:20486313..20486315,aa:Sec)"; '
        'exon_number "1";\n',
    )

    result = read_gtf(path)
    row = result.iloc[0]

    assert row["note"] == (
        "UGA stop codon recoded as selenocysteine; "
        "The RefSeq protein has 1 substitution compared to this genomic sequence"
    )
    assert row["product"] == "protein adenylyltransferase SelO, mitochondrial"
    assert row["exon_number"] == "1"
    assert row["Start"] == 20486312


def test_read_gtf_supports_unquoted_attribute_values(tmp_path: Path):
    path = _write_temp_gtf(
        tmp_path,
        "chr1\tensembl_havana\texon\t3069203\t3069296\t.\t+\t.\t"
        'gene_id "ENSG00000142611.17"; transcript_id "ENST00000270722.10"; '
        'gene_type "protein_coding"; gene_name "PRDM16"; transcript_type "protein_coding"; '
        'transcript_name "PRDM16-201"; exon_number 1; exon_id "ENSE00003850248.1"; '
        'tag "MANE_Select"; protein_id "ENSP00000270722.5"; db_xref "RefSeq:NM_022114.4";\n',
    )

    compiled = _normalize_frame_for_compare(read_gtf(path))
    python = _normalize_frame_for_compare(read_gtf_python(path))

    assert compiled.iloc[0]["exon_number"] == "1"
    assert python.iloc[0]["exon_number"] == "1"
    assert compiled.iloc[0]["protein_id"] == "ENSP00000270722.5"
    assert python.iloc[0]["protein_id"] == "ENSP00000270722.5"
    assert compiled.iloc[0]["gene_id"] == python.iloc[0]["gene_id"]
    assert compiled.iloc[0]["transcript_id"] == python.iloc[0]["transcript_id"]
    assert compiled.iloc[0]["exon_id"] == python.iloc[0]["exon_id"]


@pytest.mark.parametrize("reader", [read_gtf, read_gtf_python], ids=["default", "python"])
def test_full_readers_support_mixed_attribute_formats(tmp_path: Path, reader):
    path = _write_temp_gtf(
        tmp_path,
        "# header\n"
        "chr1\tRefSeq\tCDS\t20486313\t20486315\t.\t+\t0\t"
        'gene_id "G1"; transcript_id "T1"; note "contains; semicolon"; exon_number "1"; product "Protein 1";\n'
        "chr1\tensembl_havana\texon\t3069203\t3069296\t.\t+\t.\t"
        'gene_id "G2"; transcript_id "T2"; transcript_name "TX2"; exon_number 2; exon_id "EX2"; level 3;\n',
    )

    result = _normalize_frame_for_compare(reader(path))

    assert list(result["gene_id"]) == ["G1", "G2"]
    assert list(result["transcript_id"]) == ["T1", "T2"]
    assert result.iloc[0]["note"] == "contains; semicolon"
    assert result.iloc[0]["exon_number"] == "1"
    assert result.iloc[1]["exon_number"] == "2"
    assert result.iloc[1]["level"] == "3"
    assert result.iloc[1]["exon_id"] == "EX2"


def test_read_gtf_duplicate_attr_keeps_all_values(tmp_path: Path):
    path = _write_temp_gtf(
        tmp_path,
        "chr1\thavana\texon\t11869\t12227\t.\t+\t.\t"
        'gene_id "ENSG1"; transcript_id "ENST1"; exon_id "ENSE1"; tag "CCDS"; tag "basic";\n',
    )

    default = read_gtf(path)
    duplicate = read_gtf(path, duplicate_attr=True)

    assert default.iloc[0]["tag"] == "basic"
    assert duplicate.iloc[0]["tag"] == "CCDS,basic"


def test_read_gtf_full_multiple_chunks_match(tmp_path: Path):
    path = _write_temp_gtf(
        tmp_path,
        "# header\n"
        'chr1\thavana\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG1"; gene_name "DDX11L1"; level "2";\n'
        'chr1\thavana\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; transcript_name "DDX11L1-201"; tag "basic";\n'
        'chr1\thavana\texon\t12010\t12057\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; exon_number "1"; exon_id "ENSE1";\n',
    )

    skiprows = find_first_data_line_index(path)
    chunked = read_gtf_full(path, skiprows=skiprows, chunksize=2)
    unchunked = read_gtf_full(path, skiprows=skiprows, chunksize=1000)

    chunked = _normalize_frame_for_compare(chunked)
    unchunked = _normalize_frame_for_compare(unchunked)

    assert_frame_equal(chunked, unchunked, check_dtype=False)


def test_read_gtf_restricted_returns_core_columns(tmp_path: Path):
    path = _write_temp_gtf(
        tmp_path,
        "# header\n"
        'chr1\thavana\texon\t12010\t12057\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; exon_number "1"; exon_id "ENSE1";\n',
    )

    skiprows = find_first_data_line_index(path)
    result = read_gtf_restricted(path, skiprows=skiprows)

    assert list(result.columns) == [
        "Chromosome",
        "Source",
        "Feature",
        "Start",
        "End",
        "Score",
        "Strand",
        "Frame",
        "gene_id",
        "transcript_id",
        "exon_number",
        "exon_id",
    ]
    assert result.iloc[0]["Start"] == 12009
    assert result.iloc[0]["transcript_id"] == "ENST1"


@pytest.mark.parametrize(
    "reader",
    [read_gtf_restricted, read_gtf_restricted_python],
    ids=["default", "python"],
)
def test_restricted_readers_support_unquoted_values(tmp_path: Path, reader):
    path = _write_temp_gtf(
        tmp_path,
        "# header\n"
        'chr1\tensembl_havana\texon\t3069203\t3069296\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number 4; exon_id "EX4";\n',
    )

    skiprows = find_first_data_line_index(path)
    result = _normalize_frame_for_compare(reader(path, skiprows=skiprows))

    assert result.iloc[0]["gene_id"] == "G1"
    assert result.iloc[0]["transcript_id"] == "T1"
    assert result.iloc[0]["exon_number"] == "4"
    assert result.iloc[0]["exon_id"] == "EX4"


def test_read_gtf_python_matches_compiled_reader(tmp_path: Path):
    path = _write_temp_gtf(
        tmp_path,
        "# header\n"
        'chr1\thavana\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG1"; gene_name "DDX11L1"; level "2";\n'
        'chr1\thavana\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; transcript_name "DDX11L1-201"; tag "basic";\n',
    )

    skiprows = find_first_data_line_index(path)
    compiled = read_gtf_full(path, skiprows=skiprows, chunk_size=2)
    python = read_gtf_full_python(path, skiprows=skiprows, chunk_size=2)

    compiled = _normalize_frame_for_compare(compiled)
    python = _normalize_frame_for_compare(python)

    assert list(python["gene_id"]) == list(compiled["gene_id"])
    assert list(python["gene_name"]) == ["DDX11L1", None]
    assert list(python["transcript_id"]) == [None, "ENST1"]
    assert list(python["transcript_name"]) == [None, "DDX11L1-201"]
    assert list(python["tag"]) == [None, "basic"]
    assert list(python["Start"]) == list(compiled["Start"])


def test_read_gtf_python_duplicate_attr_keeps_all_values(tmp_path: Path):
    path = _write_temp_gtf(
        tmp_path,
        "chr1\thavana\texon\t11869\t12227\t.\t+\t.\t"
        'gene_id "ENSG1"; transcript_id "ENST1"; exon_id "ENSE1"; tag "CCDS"; tag "basic";\n',
    )

    compiled = read_gtf(path, duplicate_attr=True)
    python = read_gtf_python(path, duplicate_attr=True)

    assert compiled.iloc[0]["tag"] == "CCDS,basic"
    assert python.iloc[0]["tag"] == "CCDS,basic"


def test_read_gtf_restricted_python_matches_compiled(tmp_path: Path):
    path = _write_temp_gtf(
        tmp_path,
        "# header\n"
        'chr1\thavana\texon\t12010\t12057\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; exon_number "1"; exon_id "ENSE1";\n',
    )

    skiprows = find_first_data_line_index(path)
    compiled = read_gtf_restricted(path, skiprows=skiprows)
    python = read_gtf_restricted_python(path, skiprows=skiprows)

    compiled = _normalize_frame_for_compare(compiled)
    python = _normalize_frame_for_compare(python)

    assert_frame_equal(compiled, python, check_dtype=False)
