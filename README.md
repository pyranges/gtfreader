# gtfreader

Fast GTF reading into pandas DataFrames.

## Install

```bash
python -m pip install -e .
```

## Example

```python
from gtfreader import read_gtf

df = read_gtf("annotation.gtf")
print(df.columns)
print(df.head())
```
