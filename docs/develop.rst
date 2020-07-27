.. _develop:

Develop AvaFrame
================

Notes on design principles and how to contribute to AvaFrame


How to install
--------------

Clone repository::

  git clone https://github.com/avaframe/AvaFrame.git

Change directory::

  cd AvaFrame

And install avaframe::

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

For local documentation install sphinx with::

  pip install sphinx

(or equivalent conda command)

We use the ReadTheDocs theme, which you need to install with::

  pip install sphinx-rtd-theme

Goto docs directory and e.g.::

  make html

to generate html documentation within the _build directory.


Logging
-------

Use python provided logging module. See most simple example in runTmp1Ex.py and tmp1Ex.py

Basically: in your main script call::

  import logging
  import sys

  # create logger, set to logging.DEBUG to see all messages
  logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                    format='%(module)s:%(levelname)s - %(message)s')

And in your modules/subscripts add::

  import logging
  log = logging.getLogger(__name__)

So you can use::

  log.debug('Should be here')
  log.info('DEM is %s',variable)

To get output that looks like::

  tmp1Ex:DEBUG - Should be here 

  

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

