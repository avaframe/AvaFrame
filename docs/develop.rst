Development
===========


Notes to developers
-------------------

Here you can find notes on design principles and how to contribute to AvaFrame.

First a few general remarks. We are aware there's no right or wrong here, we
decided to stick to a few principles. If you contribute code, please try to observe
these:

AF-DESIGN-1: Use SI units
  Units in APIs should be SI units: seconds, meters, etc.

AF-STYLE-1: Line length 
  If you want a guide for line length --- we think line breaking at 120
  characters makes more sense than line breaking at 80.

AF-STYLE-2: Naming
  Use short, but descriptive enough names for variables. Use longer descriptive names for
  things with bigger scopes (functions, modules). Main theme: as short as possible, but as long 
  as necessary. 
  
  Please use lowerCamelCase (e.g.: myFancyVariable) for naming, avoid underscores.  
  
AF-STYLE-3: Commenting
  Exposed API functions need documentation comments, most importantly a docstring, 
  other code only needs to be commented as necessary. Write your code clearly and use sensible
  names to reduce the need for comments.

AF-STYLE-4: No dead code
  Unused code, commented out code, functions that are never called, etc,
  should be removed from the project to reduce the cognitive load of
  reading the source code. Old code is available in the source history if
  it is needed.

--------------------------

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

Pull requests are reviewed and handled. Once the pull request is included into the
master, the local myAwesomeFeature branch can be deleted (the one in the main
repository/origin will be handled by the pull request)::

  git checkout master
  git branch -d myAwesomeFeature


--------------------------


Build the documentation
------------------------

If you want to work on the documentation you need to install *Sphinx*. If you have followed the conda installation using
``avaframe_env_spec.txt``, you can omit the following steps. If not, you can install Sphinx, the *ReadTheDocs*
theme, and the *sphinxcontrib-bibtex*, which we use to include references, by running::

  pip install sphinx
  pip install sphinx-rtd-theme
  pip install sphinxcontrib-bibtex

In order to build the documentation you need to install make ::

  sudo apt install make

Then go to the ``docs\`` directory and run::

  make html

Html files of the documentation can be found in the _build directory.


How to test code
-----------------

AvaFrame uses pytest to test code. If you add new code, consider including a
pytest for it in ``Avaframe/avaframe/tests/``. In order to perform the pytests, just run::

  pytest

and you should see something like::

  ===========================test session starts ======================
  platform linux -- Python 3.8.3, pytest-5.4.3, py-1.9.0, pluggy-0.13.1
  collected 1 item
  tests/test_tmp1Ex.py .                                       [100%]

  ==========================1 passed in 0.02s =========================

---------------------------

How to add a benchmark test
----------------------------

AvaFrame offers an expanding benchmark test suite. At the moment this test suite
includes avalanche simulations for various idealised topographies. The
:py:mod:`runStandardTestsCom1DFA.py` facilitates running all the available benchmark tests
for com1DFA at once. With this script, the avalanche simulations are performed,
plotted and a report of the comparison between simulation results and the benchmark data is
generated. If you plan to add a new benchmark test case, follow these steps

  * first chose a name, we suggest to start it with ``ava`` (for now let's refer
    to it as ``NameOfAvalanche``)
  * add all the required input data in ``data/NameOfAvalanche``. Follow the
    required directory structure which can be generated using:
    :ref:`moduleIn3Utils:Initialize Project`

as a next step, you need to add the benchmark results:

  * go to ``AvaFrame/benchmarks`` and add the subdirectory named after your test name
  * add benchmark data i.e. peak values of result parameters as ascii files. This
    data will be used as reference for the new test!
  * add the configuration file as ``NameOfAvalanche_com1DFACfg.ini``
  * add a json file with required info on benchmark test - you can use the
    example provided in :py:mod:`runScripts/runWriteDesDict.py`
  * go to ``AvaFrame/benchmarks/simParametersDict.py`` and add a simulation
    dictionary that contains all the info on the new benchmark

Now, you are ready to go! Move to ``AvaFrame/avaframe`` and run: ::

  python runStandardTestsCom1DFA.py

You can check out the markdown-style report of the comparison at:
``tests/reports/standardTestsReportPy.md``.


