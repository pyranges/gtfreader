"""The compiled attribute parsers must match their reference implementations.

`parse_gff3_chunk_columns` must match pyranges1's `to_keys_and_values`.

That function is the reference: it is what `read_gff3` has always used, quirks
included. The reference is reimplemented here rather than imported, so gtfreader
keeps no dependency on pyranges1 and the tests still say what they are checking.
"""

from __future__ import annotations

import pytest

from gtfreader import parse_chunk_columns, parse_gff3_chunk_columns


def reference(line: str) -> dict[str, str]:
    """`pyranges1.readers.to_keys_and_values`, verbatim."""
    return dict(it.split("=", 1) for it in line.rstrip("; ").split(";") if "=" in it)


def reference_columns(lines: list[str]) -> dict[str, list]:
    """The reference applied column-wise, in first-appearance order."""
    columns: dict[str, list] = {}
    for row, line in enumerate(lines):
        for key, value in reference(line).items():
            columns.setdefault(key, [None] * len(lines))[row] = value
    return columns


CASES = [
    ("plain", ["ID=a;Name=b"]),
    ("trailing semicolon", ["ID=a;Name=b;"]),
    ("trailing semicolons and spaces", ["ID=a; ; "]),
    ("space after semicolon keeps it in the key", ["ID=a; Name=b"]),
    ("empty segment", ["ID=a;;Name=b"]),
    ("segment without an equals", ["novalue;ID=a"]),
    ("equals inside the value", ["ID=a=b=c"]),
    ("empty attribute", [""]),
    ("only semicolons", [";;;"]),
    ("repeated key keeps the last", ["ID=a;ID=b"]),
    ("empty value", ["ID=;Name=b"]),
    ("empty key", ["=v;ID=a"]),
    ("ragged rows", ["ID=a", "ID=b;Name=c", "Note=d"]),
    ("first appearance order", ["b=1;a=2", "c=3"]),
    ("utf8", ["ID=éèü;Name=ñ"]),
    ("long values", ["ID=" + "x" * 500]),
    ("many rows", [f"ID=g{i};Name=n{i % 7};tag=t{i % 3}" for i in range(1000)]),
    ("runs of identical values", ["ID=same;n=1"] * 50 + ["ID=other;n=2"] * 50),
    ("value cache overflow", [f"ID=v{i}" for i in range(20_000)]),
]


@pytest.mark.parametrize(("name", "lines"), CASES, ids=[c[0] for c in CASES])
def test_matches_the_reference(name, lines):
    assert parse_gff3_chunk_columns(lines) == reference_columns(lines)


@pytest.mark.parametrize(("name", "lines"), CASES, ids=[c[0] for c in CASES])
def test_column_order_matches_first_appearance(name, lines):
    assert list(parse_gff3_chunk_columns(lines)) == list(reference_columns(lines))


def test_every_column_is_as_long_as_the_chunk():
    lines = ["ID=a", "Name=b", "ID=c;Other=d"]
    for column in parse_gff3_chunk_columns(lines).values():
        assert len(column) == len(lines)


def test_repeated_values_are_the_same_object():
    """The point of the value cache: pandas builds a frame ~1.8x faster from
    shared objects than from equal ones."""
    columns = parse_gff3_chunk_columns(["kind=exon"] * 100 + ["kind=cds"] * 100)
    kind = columns["kind"]
    assert len({id(v) for v in kind}) == 2


def test_non_adjacent_repeats_are_also_shared():
    """The run memo alone would miss these; the dict cache catches them."""
    columns = parse_gff3_chunk_columns([f"kind={'abc'[i % 3]}" for i in range(300)])
    assert len({id(v) for v in columns["kind"]}) == 3


def test_rejects_non_strings():
    with pytest.raises(TypeError):
        parse_gff3_chunk_columns([1.5])


def test_empty_chunk():
    assert parse_gff3_chunk_columns([]) == {}


# --------------------------------------------------------------------------
# The GTF parser shares the same key/value machinery
# --------------------------------------------------------------------------


def test_gtf_repeated_values_are_the_same_object():
    columns = parse_chunk_columns(['gene_type "protein_coding";'] * 100 + ['gene_type "lncRNA";'] * 100)
    assert len({id(v) for v in columns["gene_type"]}) == 2


def test_gtf_non_adjacent_repeats_are_also_shared():
    lines = [f'tag "{"abc"[i % 3]}";' for i in range(300)]
    assert len({id(v) for v in parse_chunk_columns(lines)["tag"]}) == 3


def test_gtf_unique_values_are_not_shared():
    """A column with a distinct value per row must not fill the cache."""
    lines = [f'exon_id "ENSE{i:011d}";' for i in range(20_000)]
    assert len({id(v) for v in parse_chunk_columns(lines)["exon_id"]}) == 20_000


@pytest.mark.parametrize(
    ("name", "lines", "expected"),
    [
        ("quoted", ['gene_id "G1";'], {"gene_id": ["G1"]}),
        ("unquoted", ["exon_number 3;"], {"exon_number": ["3"]}),
        ("semicolon inside quotes", ['note "a; b";'], {"note": ["a; b"]}),
        ("repeated key keeps the last", ['tag "a"; tag "b";'], {"tag": ["b"]}),
        ("trailing whitespace trimmed", ["level 2  ;"], {"level": ["2"]}),
        ("empty", [""], {}),
        (
            "keys beyond the six that used to be special-cased",
            ['havana_gene "H1"; ccdsid "C1"; ont "O1";'],
            {"havana_gene": ["H1"], "ccdsid": ["C1"], "ont": ["O1"]},
        ),
    ],
)
def test_gtf_parses(name, lines, expected):
    assert parse_chunk_columns(lines) == expected


def test_gtf_rejects_non_strings():
    with pytest.raises(TypeError):
        parse_chunk_columns([1.5])
