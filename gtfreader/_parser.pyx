# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: initializedcheck=False
# cython: infer_types=True
# cython: nonecheck=False

"""Compiled attribute-column parsers for GTF and GFF3.

Both work on the UTF-8 bytes of each line rather than indexing the `str`.
Indexing a Python string goes through a width-aware read per character, which
measured 350-700 MB/s; the delimiters being looked for are all ASCII, so the
scan can run over the raw buffer instead and comparisons become `memcmp`.
CPython hands back the UTF-8 buffer of an ASCII string without copying, and
every GTF and GFF3 attribute in practice is ASCII; a line that is not still
works, because offsets stay byte offsets throughout and values are decoded back
with `PyUnicode_DecodeUTF8`.
"""

from libc.string cimport memcmp

cdef extern from "Python.h":
    const char* PyUnicode_AsUTF8AndSize(object unicode, Py_ssize_t *size) except NULL
    object PyUnicode_DecodeUTF8(const char *s, Py_ssize_t size, const char *errors)


cdef enum:
    TAB = 9
    SPACE = 32
    QUOTE = 34
    SEMICOLON = 59
    EQUALS = 61


# Values repeat: an annotation names the same gene on every one of its exons,
# and columns like gene_type hold a handful of values across the whole file.
# Handing pandas the *same* object each time rather than an equal one is worth
# 1.9x on frame construction at gencode-like cardinality, and frame
# construction is half of attribute expansion.
#
# Two levels, because they catch different shapes. The run memo compares
# against the column's previous value and allocates nothing at all on a hit,
# which is the common case in a coordinate-sorted file. The dict catches values
# that recur without being adjacent -- `tag` cycling through five -- but has to
# build the string before it can look it up, so it is the weaker of the two. It
# is capped: a column with a distinct value per row, like exon_id, would
# otherwise fill it with a million entries it can never hit.
DEF VALUE_CACHE_CAP = 8192


cdef inline bint _matches(const char* s, Py_ssize_t a, Py_ssize_t b, object other):
    """Is s[a:b] equal to `other`, without building s[a:b]?"""
    cdef Py_ssize_t n = b - a
    cdef Py_ssize_t other_len
    cdef const char* other_buf

    if other is None:
        return False
    other_buf = PyUnicode_AsUTF8AndSize(other, &other_len)
    if other_len != n:
        return False
    if n == 0:
        return True
    return memcmp(s + a, other_buf, n) == 0


cdef inline Py_ssize_t _key_index(
    const char* s,
    Py_ssize_t a,
    Py_ssize_t b,
    list keys,
    Py_ssize_t ordinal,
):
    """Index of s[a:b] in `keys`, or -1. Does not build s[a:b].

    A file has a handful of distinct attribute keys, repeats them on every row,
    and -- this is the part worth exploiting -- lists them in the same order
    every time. So the key in position `ordinal` of this row is almost always
    the one that was in position `ordinal` of the last row: try that first and
    the scan is one comparison, not half the table.
    """
    cdef Py_ssize_t i
    cdef Py_ssize_t n = len(keys)

    if 0 <= ordinal < n and _matches(s, a, b, <str> keys[ordinal]):
        return ordinal
    for i in range(n):
        if _matches(s, a, b, <str> keys[i]):
            return i
    return -1


cdef inline object _value_at(
    const char* s,
    Py_ssize_t a,
    Py_ssize_t b,
    list last_values,
    list caches,
    Py_ssize_t index,
):
    """s[a:b], reusing an equal object already held for this column."""
    cdef object value, cached, last
    cdef dict cache

    last = last_values[index]
    if _matches(s, a, b, last):
        return last

    value = PyUnicode_DecodeUTF8(s + a, b - a, NULL)
    cache = <dict> caches[index]
    cached = cache.get(value)
    if cached is None:
        if len(cache) < VALUE_CACHE_CAP:
            cache[value] = value
    else:
        value = cached

    last_values[index] = value
    return value


cdef inline Py_ssize_t _add_column(
    object key,
    Py_ssize_t rows,
    dict columns,
    list keys,
    list cols,
    list last_values,
    list caches,
):
    cdef list col = [None] * rows

    keys.append(key)
    cols.append(col)
    last_values.append(None)
    caches.append({})
    columns[key] = col
    return len(keys) - 1


def parse_chunk_columns(lines):
    """Expand GTF `key "value";` attribute strings into a dict of columns.

    Column-oriented: one list per key, filled in place, so the caller can hand
    the dict straight to pandas without building a dict per row first.
    """
    cdef Py_ssize_t n, row_i, m
    cdef Py_ssize_t pos, key_start, key_end, val_start, val_end, index, ordinal
    cdef dict columns
    cdef list keys, cols, last_values, caches, col
    cdef str line
    cdef const char* s

    n = len(lines)
    columns = {}
    keys = []
    cols = []
    last_values = []
    caches = []

    for row_i in range(n):
        line = lines[row_i]
        s = PyUnicode_AsUTF8AndSize(line, &m)
        pos = 0
        ordinal = 0

        while pos < m:
            while pos < m and (s[pos] == SPACE or s[pos] == TAB or s[pos] == SEMICOLON):
                pos += 1
            if pos >= m:
                break

            key_start = pos

            while pos < m and s[pos] != SPACE and s[pos] != TAB:
                pos += 1
            key_end = pos

            while pos < m and (s[pos] == SPACE or s[pos] == TAB):
                pos += 1
            if pos >= m:
                break

            if s[pos] == QUOTE:
                val_start = pos + 1
                pos += 1

                while pos < m and s[pos] != QUOTE:
                    pos += 1
                if pos >= m:
                    break

                val_end = pos
            else:
                val_start = pos
                while pos < m and s[pos] != SEMICOLON:
                    pos += 1
                val_end = pos
                while val_end > val_start and (s[val_end - 1] == SPACE or s[val_end - 1] == TAB):
                    val_end -= 1

            index = _key_index(s, key_start, key_end, keys, ordinal)
            if index < 0:
                index = _add_column(
                    PyUnicode_DecodeUTF8(s + key_start, key_end - key_start, NULL),
                    n,
                    columns,
                    keys,
                    cols,
                    last_values,
                    caches,
                )
            col = <list> cols[index]

            col[row_i] = _value_at(s, val_start, val_end, last_values, caches, index)

            ordinal += 1
            pos += 1

    return columns


def parse_gff3_chunk_columns(lines):
    """Expand GFF3 `key=value;` attribute strings into a dict of columns.

    Matches `pyranges1.readers.to_keys_and_values` exactly, quirks included:

      - the string is right-stripped of `;` and spaces, then split on `;`
      - a segment with no `=` is skipped, which is also what makes an empty
        attribute an empty row rather than an error
      - the split is at the *first* `=`, so a value may contain more
      - no whitespace is trimmed around a key, so `ID=a; Name=b` really does
        yield a key of `" Name"`. GFF3 does not put a space there, but files do
      - a repeated key keeps the last value
      - columns come out in order of first appearance across the chunk

    Column-oriented on purpose. Building a dict per row and handing the lot to
    `DataFrame.from_records` costs as much again as the parsing does; filling
    one list per column skips that entirely.
    """
    cdef Py_ssize_t n, row_i, m
    cdef Py_ssize_t pos, end, seg_start, seg_end, eq, index, ordinal
    cdef dict columns
    cdef list keys, cols, last_values, caches, col
    cdef str line
    cdef const char* s

    n = len(lines)
    columns = {}
    keys = []
    cols = []
    last_values = []
    caches = []

    for row_i in range(n):
        line = lines[row_i]
        s = PyUnicode_AsUTF8AndSize(line, &m)

        # `line.rstrip("; ")`
        end = m
        while end > 0 and (s[end - 1] == SEMICOLON or s[end - 1] == SPACE):
            end -= 1

        pos = 0
        ordinal = 0
        while pos < end:
            seg_start = pos
            while pos < end and s[pos] != SEMICOLON:
                pos += 1
            seg_end = pos
            pos += 1

            eq = seg_start
            while eq < seg_end and s[eq] != EQUALS:
                eq += 1
            if eq >= seg_end:
                # No `=` in this segment: not a tag=value pair.
                continue

            index = _key_index(s, seg_start, eq, keys, ordinal)
            if index < 0:
                index = _add_column(
                    PyUnicode_DecodeUTF8(s + seg_start, eq - seg_start, NULL),
                    n,
                    columns,
                    keys,
                    cols,
                    last_values,
                    caches,
                )
            col = <list> cols[index]

            col[row_i] = _value_at(s, eq + 1, seg_end, last_values, caches, index)
            ordinal += 1

    return columns
