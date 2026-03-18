# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: initializedcheck=False
# cython: infer_types=True
# cython: nonecheck=False

cdef inline int _key_code(str s, Py_ssize_t a, Py_ssize_t b):
    cdef Py_ssize_t n
    n = b - a

    if n == 7:
        if (
            s[a] == 'g' and s[a + 1] == 'e' and s[a + 2] == 'n' and
            s[a + 3] == 'e' and s[a + 4] == '_' and s[a + 5] == 'i' and
            s[a + 6] == 'd'
        ):
            return 1
        if (
            s[a] == 'e' and s[a + 1] == 'x' and s[a + 2] == 'o' and
            s[a + 3] == 'n' and s[a + 4] == '_' and s[a + 5] == 'i' and
            s[a + 6] == 'd'
        ):
            return 6

    elif n == 9:
        if (
            s[a] == 'g' and s[a + 1] == 'e' and s[a + 2] == 'n' and
            s[a + 3] == 'e' and s[a + 4] == '_' and s[a + 5] == 'n' and
            s[a + 6] == 'a' and s[a + 7] == 'm' and s[a + 8] == 'e'
        ):
            return 3

    elif n == 11:
        if (
            s[a] == 'e' and s[a + 1] == 'x' and s[a + 2] == 'o' and
            s[a + 3] == 'n' and s[a + 4] == '_' and s[a + 5] == 'n' and
            s[a + 6] == 'u' and s[a + 7] == 'm' and s[a + 8] == 'b' and
            s[a + 9] == 'e' and s[a + 10] == 'r'
        ):
            return 5

    elif n == 13:
        if (
            s[a] == 't' and s[a + 1] == 'r' and s[a + 2] == 'a' and
            s[a + 3] == 'n' and s[a + 4] == 's' and s[a + 5] == 'c' and
            s[a + 6] == 'r' and s[a + 7] == 'i' and s[a + 8] == 'p' and
            s[a + 9] == 't' and s[a + 10] == '_' and s[a + 11] == 'i' and
            s[a + 12] == 'd'
        ):
            return 2

    elif n == 15:
        if (
            s[a] == 't' and s[a + 1] == 'r' and s[a + 2] == 'a' and
            s[a + 3] == 'n' and s[a + 4] == 's' and s[a + 5] == 'c' and
            s[a + 6] == 'r' and s[a + 7] == 'i' and s[a + 8] == 'p' and
            s[a + 9] == 't' and s[a + 10] == '_' and s[a + 11] == 'n' and
            s[a + 12] == 'a' and s[a + 13] == 'm' and s[a + 14] == 'e'
        ):
            return 4

    return 0


cdef inline str _cached_str(dict cache, str value):
    cdef object obj
    obj = cache.get(value)
    if obj is None:
        cache[value] = value
        return value
    return <str> obj


def parse_chunk_columns(lines):
    cdef Py_ssize_t n, row_i
    cdef Py_ssize_t pos, m, key_start, key_end, val_start, val_end
    cdef dict columns
    cdef dict gene_name_cache

    cdef str s, key, value
    cdef list col
    cdef int code

    cdef object obj_col

    cdef list gene_id_col
    cdef list transcript_id_col
    cdef list gene_name_col
    cdef list transcript_name_col
    cdef list exon_number_col
    cdef list exon_id_col

    n = len(lines)
    columns = {}
    gene_name_cache = {}

    gene_id_col = None
    transcript_id_col = None
    gene_name_col = None
    transcript_name_col = None
    exon_number_col = None
    exon_id_col = None

    for row_i in range(n):
        s = <str> lines[row_i]
        m = len(s)
        pos = 0

        while pos < m:
            while pos < m and (s[pos] == ' ' or s[pos] == '\t' or s[pos] == ';'):
                pos += 1
            if pos >= m:
                break

            key_start = pos

            while pos < m and s[pos] != ' ' and s[pos] != '\t':
                pos += 1
            key_end = pos

            while pos < m and (s[pos] == ' ' or s[pos] == '\t'):
                pos += 1
            if pos >= m:
                break

            if s[pos] == '"':
                val_start = pos + 1
                pos += 1

                while pos < m and s[pos] != '"':
                    pos += 1
                if pos >= m:
                    break

                val_end = pos
            else:
                val_start = pos
                while pos < m and s[pos] != ';':
                    pos += 1
                val_end = pos
                while val_end > val_start and (s[val_end - 1] == ' ' or s[val_end - 1] == '\t'):
                    val_end -= 1

            value = s[val_start:val_end]

            code = _key_code(s, key_start, key_end)

            if code == 1:
                if gene_id_col is None:
                    gene_id_col = [None] * n
                    columns["gene_id"] = gene_id_col
                gene_id_col[row_i] = value

            elif code == 2:
                if transcript_id_col is None:
                    transcript_id_col = [None] * n
                    columns["transcript_id"] = transcript_id_col
                transcript_id_col[row_i] = value

            elif code == 3:
                if gene_name_col is None:
                    gene_name_col = [None] * n
                    columns["gene_name"] = gene_name_col
                gene_name_col[row_i] = _cached_str(gene_name_cache, value)

            elif code == 4:
                if transcript_name_col is None:
                    transcript_name_col = [None] * n
                    columns["transcript_name"] = transcript_name_col
                transcript_name_col[row_i] = value

            elif code == 5:
                if exon_number_col is None:
                    exon_number_col = [None] * n
                    columns["exon_number"] = exon_number_col
                exon_number_col[row_i] = value

            elif code == 6:
                if exon_id_col is None:
                    exon_id_col = [None] * n
                    columns["exon_id"] = exon_id_col
                exon_id_col[row_i] = value

            else:
                key = s[key_start:key_end]
                obj_col = columns.get(key)
                if obj_col is None:
                    col = [None] * n
                    columns[key] = col
                else:
                    col = <list> obj_col
                col[row_i] = value

            pos += 1

    return columns
