Develop AvaFrame
================

Here you can find notes on design principles and how to contribute to AvaFrame.


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
pytest for it in `Avaframe/avaframe/tests/``. In order to perform the pytests, just run::

  pytest

and you should see something like::

  ===========================test session starts ======================
  platform linux -- Python 3.8.3, pytest-5.4.3, py-1.9.0, pluggy-0.13.1
  collected 1 item
  tests/test_tmp1Ex.py .                                       [100%]

  ==========================1 passed in 0.02s =========================
