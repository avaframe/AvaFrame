# AvaFrame

This is the main AvaFrame repository

## Directory structure:

```
-avaframe: main python scripts etc
    |-data: data needed for calibration/tests/etc
    |-tests: pytest scripts

-benchmarks: our references 

-docs: md/rst for ReadTheDoc documentation
```


## Documentation

For the documentation install sphinx with 

```
pip install sphinx
```
(or equivalent conda command)

We use the ReadTheDocs theme, which you need to install with

```
pip install sphinx-rtd-theme
```

Goto docs directory and e.g.

```
make html
```

to generate html documentation within the _build directory
