# gtfread

`gtfread` is a small Cython-backed package for parsing GTF attribute fields into column-oriented Python lists.

## Layout

- `src/gtfread/_parser.pyx`: the Cython implementation
- `src/gtfread/_fallback.py`: pure Python fallback for source-tree imports before the extension is built
- `tests/test_parser.py`: smoke test for the public API

## Install

```bash
python -m pip install -e .
```

That uses `pyproject.toml` to install the build requirements and compile the Cython extension.

## Build a wheel

```bash
python -m build
```

## In-place extension build

```bash
python setup.py build_ext --inplace
```

## Do you need `setup.py`?

No. A modern project can build a Cython extension from `pyproject.toml` alone if you declare the extension there and use a build backend that supports it. `setup.py` is still useful when you want `cythonize(...)` directly, custom build logic, or more control over compiler options.
