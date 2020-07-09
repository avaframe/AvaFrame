# This is the AvaFrame repository

# Installation, documentation

The documentation is hosted on ReadTheDocs: http://docs.avaframe.org

For local documentation install sphinx with:

```
  pip install sphinx
```

(or equivalent conda command)

We use the ReadTheDocs theme, which you need to install with::

```
  pip install sphinx-rtd-theme
```

Goto docs directory and e.g.::

```
  make html
```

to generate html documentation within the _build directory.




# Directory structure:

```
  -avaframe: main python scripts etc
      |-data: data needed for calibration/tests/etc
      |-tests: pytest scripts

  -benchmarks: our references 

  -docs: rst for ReadTheDoc documentation
```



# Tests 

[<img src="https://readthedocs.org/projects/avaframe/badge/?version=latest">](http://docs.avaframe.org/en/latest/)

# License 
Licensed with [![European Public License EUPL](https://img.shields.io/badge/license-EUPL-green.png)](https://git.avaframe.org/AvaFrame/AvaFrame/src/branch/master/LICENSE.txt)

