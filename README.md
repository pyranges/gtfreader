# gtfreader

Fast GTF reading into pandas DataFrames.

## Install

```bash
python -m pip install -e .
```

### Optional: parse on every core

pandas' CSV parser is single-threaded, and on a large GTF the nine fixed
columns are about a third of the cost of reading. With `pyarrow` installed that
part is parsed on every core instead:

```bash
python -m pip install -e ".[fast-io]"
```

`read_gtf` is 1.5x faster end to end at 10^6 rows on twelve cores. It is only
the tabular parse that speeds up; expanding the attribute column is the rest of
the work and is unchanged.

pandas remains the reference implementation and the fallback. The pyarrow path
is skipped when pyarrow is missing, when `nrows` is set, when a file is shaped
in a way it does not model, and whenever the two parsers would disagree -- the
result is the same frame either way.

## Example

```python
from gtfreader import read_gtf

df = read_gtf("annotation.gtf")
print(df.columns)
print(df.head())
```
