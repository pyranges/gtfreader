# gtfreader

`gtfreader` is a small package for parsing and reading GTF files into pandas dataframes.

## Install

```bash
python -m pip install -e .
```

## Usage

```python
from gtfreader import read_gtf, read_gtf_python

df = read_gtf("annotation.gtf")
df_python = read_gtf_python("annotation.gtf")
```

`read_gtf(...)` uses the compiled parser path when available. `read_gtf_python(...)` uses the high-level pure Python parser path used by the current `pyrunges` reader style.

If you want to use the compiled low-level parser directly, pass it raw attribute strings from column 9 of the GTF before they have been expanded:

```python
import pandas as pd

from gtfreader import find_first_data_line_index, parse_chunk_columns

skiprows = find_first_data_line_index("annotation.gtf")
attribute_lines = pd.read_csv(
    "annotation.gtf",
    sep="\t",
    header=None,
    usecols=[8],
    names=["Attribute"],
    comment="#",
    skiprows=skiprows,
)["Attribute"].tolist()

compiled_columns = parse_chunk_columns(attribute_lines)
```

## Build

```bash
python -m build
```

## Test

```bash
python -m pytest -q
```
