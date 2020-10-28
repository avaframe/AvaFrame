Develop AvaFrame
================

Notes on design principles how to use and how to contribute to AvaFrame


How to install
--------------

Clone repository::

  git clone https://github.com/avaframe/AvaFrame.git


Avaframe relies on the following python library:

* numpy
* matplotlib
* pyshp
* scipy
* cmocean
* seaborn
* pytest

We suggest to work with conda and create a specific avaframe
`environment <https://conda.io/projects/conda/en/latest/user-guide/concepts/environments.html>`_ ::

    conda create --name avaframe_env

Don't forget to activate it before going on::

    conda activate avaframe_env

Move to ``Avaframe/`` repository and install the different modules::

    conda install --name avaframe_env --file avaframe_env_spec.txt

Or the equivalent pip commands.

Finally install avaframe::

  conda develop .

Or if you are using pip ::

  pip install -e .

This installs avaframe in editable mode, so everytime you import avaframe the
current (local) version will be used.


How to test
-----------

To run some simple tests just run::

  pytest

and you should see something like::

  ===========================test session starts ======================
  platform linux -- Python 3.8.3, pytest-5.4.3, py-1.9.0, pluggy-0.13.1
  collected 1 item
  tests/test_tmp1Ex.py .                                       [100%]

  ==========================1 passed in 0.02s =========================


Local documentation
-------------------

For local documentation install sphinx (already included if you used the
conda commands) with::

  pip install sphinx

We use the ReadTheDocs theme (also included in the conda environment),
which you need to install with::

  pip install sphinx-rtd-theme

We use ``sphinxcontrib-bibtex`` for references. Install it with ::

  pip install sphinxcontrib-bibtex

Go to now go to ``docs\`` directory and e.g.::

  make html

to generate html documentation within the _build directory.


Logging
-------

.. include:: moduleLogging.rst


Configuration
--------------

.. include:: moduleConfiguration.rst

Our suggested git workflow
--------------------------

Clone repository::

  git clone https://github.com/avaframe/AvaFrame.git

Clones the repository to your local machine into the directory AvaFrame. Sets
the repository to track to *origin*

Branch::

  git checkout -b myAwesomeFeature

This changes your working directory to the myAwesomeFeature branch. Try to keep
any changes in this branch specific to one bug or feature. You can have many
branches and switch in between them using the git checkout command.

Work on it and from time to time commit your changes using following commands as
necessary::

  git add
  git commit

To update this branch, you need to retrieve the changes from the master branch::

  git rebase origin master

or::

  git checkout master
  git pull
  git checkout myAwesomeFeature
  git rebase master

This replays all your changes on the current status of the master (i.e main)
branch. If conflicts arise, now is the time to solve them.

Push your changes to the main repository::

  git push origin

Once you feel you are done, start a pull request on  github.com_.

.. _github.com: https://github.com/avaframe/AvaFrame

Pull request are reviewed and handled. Once the pull request is included into the
master, the local myAwesomeFeature branch can be deleted (the one in the main
repository/origin will be handled by the pull request)::

  git checkout master
  git branch -d myAwesomeFeature
