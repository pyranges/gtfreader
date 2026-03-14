"""Pure Python fallback used when the compiled extension is not available."""


def _key_code(s: str, a: int, b: int) -> int:
    n = b - a

    if n == 3:
        if s[a] == "t" and s[a + 1] == "a" and s[a + 2] == "g":
            return 9
        if s[a] == "o" and s[a + 1] == "n" and s[a + 2] == "t":
            return 17

    elif n == 5:
        if (
            s[a] == "l"
            and s[a + 1] == "e"
            and s[a + 2] == "v"
            and s[a + 3] == "e"
            and s[a + 4] == "l"
        ):
            return 16

    elif n == 6:
        if (
            s[a] == "c"
            and s[a + 1] == "c"
            and s[a + 2] == "d"
            and s[a + 3] == "s"
            and s[a + 4] == "i"
            and s[a + 5] == "d"
        ):
            return 14

    elif n == 7:
        if (
            s[a] == "g"
            and s[a + 1] == "e"
            and s[a + 2] == "n"
            and s[a + 3] == "e"
            and s[a + 4] == "_"
            and s[a + 5] == "i"
            and s[a + 6] == "d"
        ):
            return 1
        if (
            s[a] == "e"
            and s[a + 1] == "x"
            and s[a + 2] == "o"
            and s[a + 3] == "n"
            and s[a + 4] == "_"
            and s[a + 5] == "i"
            and s[a + 6] == "d"
        ):
            return 8
        if (
            s[a] == "h"
            and s[a + 1] == "g"
            and s[a + 2] == "n"
            and s[a + 3] == "c"
            and s[a + 4] == "_"
            and s[a + 5] == "i"
            and s[a + 6] == "d"
        ):
            return 13

    elif n == 9:
        if (
            s[a] == "g"
            and s[a + 1] == "e"
            and s[a + 2] == "n"
            and s[a + 3] == "e"
            and s[a + 4] == "_"
        ):
            if (
                s[a + 5] == "n"
                and s[a + 6] == "a"
                and s[a + 7] == "m"
                and s[a + 8] == "e"
            ):
                return 3
            if (
                s[a + 5] == "t"
                and s[a + 6] == "y"
                and s[a + 7] == "p"
                and s[a + 8] == "e"
            ):
                return 4

    elif n == 10:
        if (
            s[a] == "a"
            and s[a + 1] == "r"
            and s[a + 2] == "t"
            and s[a + 3] == "i"
            and s[a + 4] == "f"
            and s[a + 5] == "_"
            and s[a + 6] == "d"
            and s[a + 7] == "u"
            and s[a + 8] == "p"
            and s[a + 9] == "l"
        ):
            return 15

    elif n == 11:
        if (
            s[a] == "e"
            and s[a + 1] == "x"
            and s[a + 2] == "o"
            and s[a + 3] == "n"
            and s[a + 4] == "_"
            and s[a + 5] == "n"
            and s[a + 6] == "u"
            and s[a + 7] == "m"
            and s[a + 8] == "b"
            and s[a + 9] == "e"
            and s[a + 10] == "r"
        ):
            return 7
        if (
            s[a] == "h"
            and s[a + 1] == "a"
            and s[a + 2] == "v"
            and s[a + 3] == "a"
            and s[a + 4] == "n"
            and s[a + 5] == "a"
            and s[a + 6] == "_"
            and s[a + 7] == "g"
            and s[a + 8] == "e"
            and s[a + 9] == "n"
            and s[a + 10] == "e"
        ):
            return 11

    elif n == 13:
        if (
            s[a] == "t"
            and s[a + 1] == "r"
            and s[a + 2] == "a"
            and s[a + 3] == "n"
            and s[a + 4] == "s"
            and s[a + 5] == "c"
            and s[a + 6] == "r"
            and s[a + 7] == "i"
            and s[a + 8] == "p"
            and s[a + 9] == "t"
            and s[a + 10] == "_"
            and s[a + 11] == "i"
            and s[a + 12] == "d"
        ):
            return 2

    elif n == 15:
        if (
            s[a] == "t"
            and s[a + 1] == "r"
            and s[a + 2] == "a"
            and s[a + 3] == "n"
            and s[a + 4] == "s"
            and s[a + 5] == "c"
            and s[a + 6] == "r"
            and s[a + 7] == "i"
            and s[a + 8] == "p"
            and s[a + 9] == "t"
            and s[a + 10] == "_"
        ):
            if (
                s[a + 11] == "n"
                and s[a + 12] == "a"
                and s[a + 13] == "m"
                and s[a + 14] == "e"
            ):
                return 5
            if (
                s[a + 11] == "t"
                and s[a + 12] == "y"
                and s[a + 13] == "p"
                and s[a + 14] == "e"
            ):
                return 6

    elif n == 17:
        if (
            s[a] == "h"
            and s[a + 1] == "a"
            and s[a + 2] == "v"
            and s[a + 3] == "a"
            and s[a + 4] == "n"
            and s[a + 5] == "a"
            and s[a + 6] == "_"
            and s[a + 7] == "t"
            and s[a + 8] == "r"
            and s[a + 9] == "a"
            and s[a + 10] == "n"
            and s[a + 11] == "s"
            and s[a + 12] == "c"
            and s[a + 13] == "r"
            and s[a + 14] == "i"
            and s[a + 15] == "p"
            and s[a + 16] == "t"
        ):
            return 10

    elif n == 24:
        if (
            s[a] == "t"
            and s[a + 1] == "r"
            and s[a + 2] == "a"
            and s[a + 3] == "n"
            and s[a + 4] == "s"
            and s[a + 5] == "c"
            and s[a + 6] == "r"
            and s[a + 7] == "i"
            and s[a + 8] == "p"
            and s[a + 9] == "t"
            and s[a + 10] == "_"
            and s[a + 11] == "s"
            and s[a + 12] == "u"
            and s[a + 13] == "p"
            and s[a + 14] == "p"
            and s[a + 15] == "o"
            and s[a + 16] == "r"
            and s[a + 17] == "t"
            and s[a + 18] == "_"
            and s[a + 19] == "l"
            and s[a + 20] == "e"
            and s[a + 21] == "v"
            and s[a + 22] == "e"
            and s[a + 23] == "l"
        ):
            return 12

    return 0


def _cached_str(cache: dict[str, str], value: str) -> str:
    cached_value = cache.get(value)
    if cached_value is None:
        cache[value] = value
        return value
    return cached_value


def parse_chunk_columns(lines):
    n = len(lines)
    columns = {}

    gene_name_cache = {}
    gene_type_cache = {}
    transcript_type_cache = {}
    tag_cache = {}
    level_cache = {}
    ont_cache = {}
    transcript_support_level_cache = {}
    artif_dupl_cache = {}

    gene_id_col = [None] * n
    transcript_id_col = [None] * n
    gene_name_col = [None] * n
    gene_type_col = [None] * n
    transcript_name_col = [None] * n
    transcript_type_col = [None] * n
    exon_number_col = [None] * n
    exon_id_col = [None] * n
    tag_col = [None] * n
    havana_transcript_col = [None] * n
    havana_gene_col = [None] * n
    transcript_support_level_col = [None] * n
    hgnc_id_col = [None] * n
    ccdsid_col = [None] * n
    artif_dupl_col = [None] * n
    level_col = [None] * n
    ont_col = [None] * n

    columns["gene_id"] = gene_id_col
    columns["transcript_id"] = transcript_id_col
    columns["gene_name"] = gene_name_col
    columns["gene_type"] = gene_type_col
    columns["transcript_name"] = transcript_name_col
    columns["transcript_type"] = transcript_type_col
    columns["exon_number"] = exon_number_col
    columns["exon_id"] = exon_id_col
    columns["tag"] = tag_col
    columns["havana_transcript"] = havana_transcript_col
    columns["havana_gene"] = havana_gene_col
    columns["transcript_support_level"] = transcript_support_level_col
    columns["hgnc_id"] = hgnc_id_col
    columns["ccdsid"] = ccdsid_col
    columns["artif_dupl"] = artif_dupl_col
    columns["level"] = level_col
    columns["ont"] = ont_col

    for row_i, s in enumerate(lines):
        m = len(s)
        pos = 0

        while pos < m:
            while pos < m and (s[pos] == " " or s[pos] == "\t" or s[pos] == ";"):
                pos += 1
            if pos >= m:
                break

            key_start = pos
            while pos < m and s[pos] != " " and s[pos] != "\t":
                pos += 1
            key_end = pos

            while pos < m and s[pos] != '"':
                pos += 1
            if pos >= m:
                break

            val_start = pos + 1
            pos += 1

            while pos < m and s[pos] != '"':
                pos += 1
            if pos >= m:
                break

            val_end = pos
            value = s[val_start:val_end]
            code = _key_code(s, key_start, key_end)

            if code == 1:
                gene_id_col[row_i] = value
            elif code == 2:
                transcript_id_col[row_i] = value
            elif code == 3:
                gene_name_col[row_i] = _cached_str(gene_name_cache, value)
            elif code == 4:
                gene_type_col[row_i] = _cached_str(gene_type_cache, value)
            elif code == 5:
                transcript_name_col[row_i] = value
            elif code == 6:
                transcript_type_col[row_i] = _cached_str(transcript_type_cache, value)
            elif code == 7:
                exon_number_col[row_i] = value
            elif code == 8:
                exon_id_col[row_i] = value
            elif code == 9:
                tag_col[row_i] = _cached_str(tag_cache, value)
            elif code == 10:
                havana_transcript_col[row_i] = value
            elif code == 11:
                havana_gene_col[row_i] = value
            elif code == 12:
                transcript_support_level_col[row_i] = _cached_str(
                    transcript_support_level_cache, value
                )
            elif code == 13:
                hgnc_id_col[row_i] = value
            elif code == 14:
                ccdsid_col[row_i] = value
            elif code == 15:
                artif_dupl_col[row_i] = _cached_str(artif_dupl_cache, value)
            elif code == 16:
                level_col[row_i] = _cached_str(level_cache, value)
            elif code == 17:
                ont_col[row_i] = _cached_str(ont_cache, value)
            else:
                key = s[key_start:key_end]
                col = columns.get(key)
                if col is None:
                    col = [None] * n
                    columns[key] = col
                col[row_i] = value

            pos += 1

    return columns
