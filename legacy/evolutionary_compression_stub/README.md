# evolutionary-compression (legacy shim)

This package exists solely to help users migrate from the historical
`evolutionary-compression` distribution on PyPI to the new name, `ecomp`.

Installing `evolutionary-compression` now emits a deprecation warning and
re-exports the `ecomp` public API. New projects should depend on `ecomp`
explicitly:

```bash
pip install ecomp
```
