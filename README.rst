
This is the AvaFrame repository
-------------------------------

Installation, documentation
---------------------------

The documentation is hosted on ReadTheDocs: http://docs.avaframe.org

For local documentation install sphinx with::

  pip install sphinx

(or equivalent conda command)

We use the ReadTheDocs theme, which you need to install with::

  pip install sphinx-rtd-theme

Goto docs directory and e.g.::

  make html

to generate html documentation within the _build directory.




Directory structure:
--------------------

::

  -avaframe: main python scripts etc
      |-data: data needed for calibration/tests/etc
      |-tests: pytest scripts

  -benchmarks: our references 

  -docs: rst for ReadTheDoc documentation





About
-----

.. :Citation:
..      .. image:: https://zenodo.org/badge/43965645.svg
..        :target: https://zenodo.org/badge/latestdoi/43965645
..        :alt: Zenodo

:Tests:       
    .. image:: https://readthedocs.org/projects/avaframe/badge/?version=latest
        :target: http://docs.avaframe.org/en/latest/
        :alt: Documentation status

..    .. image:: https://img.shields.io/badge/benchmarked%20by-asv-green.svg?style=flat

:License:
    .. image:: https://img.shields.io/badge/license%20EUPL-green.svg?style=flat
        :target: https://git.avaframe.org/AvaFrame/AvaFrame/src/branch/master/LICENSE.txt
        :alt: European Union Public License (EUPL) 

.. :Authors:

..    See the `version history`_ for a list of all contributors.

..   .. _version history: http://docs.oggm.org/en/latest/whats-new.html
