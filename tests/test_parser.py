import unittest

from gtfread import parse_chunk_columns


class ParseChunkColumnsTests(unittest.TestCase):
    def test_known_and_dynamic_columns_are_parsed(self):
        lines = [
            (
                'gene_id "GENE1"; transcript_id "TX1"; gene_name "ABC1"; '
                'gene_type "protein_coding"; transcript_name "ABC1-201"; '
                'transcript_type "protein_coding"; exon_number "1"; exon_id "EXON1"; '
                'tag "basic"; havana_transcript "OTTT0001"; havana_gene "OTTG0001"; '
                'transcript_support_level "1"; hgnc_id "HGNC:5"; ccdsid "CCDS1"; '
                'artif_dupl "false"; level "2"; ont "PGO:0000001"; custom_key "custom";'
            ),
            'gene_id "GENE2"; gene_name "ABC1"; gene_type "protein_coding"; level "2";',
        ]

        columns = parse_chunk_columns(lines)

        self.assertEqual(columns["gene_id"], ["GENE1", "GENE2"])
        self.assertEqual(columns["transcript_id"], ["TX1", None])
        self.assertEqual(columns["havana_gene"], ["OTTG0001", None])
        self.assertEqual(columns["hgnc_id"], ["HGNC:5", None])
        self.assertEqual(columns["custom_key"], ["custom", None])
        self.assertIs(columns["gene_name"][0], columns["gene_name"][1])
        self.assertIs(columns["gene_type"][0], columns["gene_type"][1])
        self.assertIs(columns["level"][0], columns["level"][1])


if __name__ == "__main__":
    unittest.main()
