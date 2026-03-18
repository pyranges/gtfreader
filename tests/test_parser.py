import pytest

parse_chunk_columns = pytest.importorskip("gtfreader._parser").parse_chunk_columns


def test_known_and_dynamic_columns_are_parsed():
    lines = [
        (
            'gene_id "GENE1"; transcript_id "TX1"; gene_name "ABC1"; '
            'gene_type "protein_coding"; transcript_name "ABC1-201"; '
            'transcript_type "protein_coding"; exon_number 1; exon_id "EXON1"; '
            'tag "basic"; havana_transcript "OTTT0001"; havana_gene "OTTG0001"; '
            'transcript_support_level "1"; hgnc_id "HGNC:5"; ccdsid "CCDS1"; '
            'artif_dupl "false"; level "2"; ont "PGO:0000001"; custom_key "custom";'
        ),
        'gene_id "GENE2"; gene_name "ABC1"; gene_type "protein_coding"; level "2";',
    ]

    columns = parse_chunk_columns(lines)

    assert columns["gene_id"] == ["GENE1", "GENE2"]
    assert columns["transcript_id"] == ["TX1", None]
    assert columns["havana_gene"] == ["OTTG0001", None]
    assert columns["hgnc_id"] == ["HGNC:5", None]
    assert columns["exon_number"] == ["1", None]
    assert columns["custom_key"] == ["custom", None]
    assert columns["gene_name"][0] is columns["gene_name"][1]
    assert columns["gene_type"] == ["protein_coding", "protein_coding"]
    assert columns["level"] == ["2", "2"]


def test_compiled_parser_omits_unseen_fast_path_columns():
    columns = parse_chunk_columns(['gene_id "GENE1"; transcript_id "TX1";'])

    assert set(columns) == {"gene_id", "transcript_id"}


def test_compiled_parser_supports_semicolons_in_quotes_and_unquoted_values():
    lines = [
        'gene_id "GENE1"; transcript_id "TX1"; note "contains; semicolon"; exon_number 2; custom_key custom_value;',
    ]

    columns = parse_chunk_columns(lines)

    assert columns["gene_id"] == ["GENE1"]
    assert columns["transcript_id"] == ["TX1"]
    assert columns["note"] == ["contains; semicolon"]
    assert columns["exon_number"] == ["2"]
    assert columns["custom_key"] == ["custom_value"]
