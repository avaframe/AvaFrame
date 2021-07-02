Installation and first run
==========================

.. Note::
  There are currently two ways to install and use AvaFrame:
   * If you want to use the operational workflow, running *com1DFA* (dense flow
     avalanche) and *com2AB* (alpha beta) with standard settings from within
     *QGis*, head over to the :ref:`installation:Operational setup and run`
     section. Use this if you are on *Windows*.
   * If you want to build your own workflows and experiment with all modules,
     head over to the :ref:`installation:Full setup and run` section. Use this
     if you are familiar with programming in python and the terminal.

.. include:: operational.rst

-----

Full setup and run
------------------

This is a quick guide on how to install AvaFrame
and the required dependencies on your machine. AvaFrame is developed on **Linux
machines** (Ubuntu/Manjaro/Arch) with recent Python versions 3.8/3.9.
These instructions assume you are familiar with working in a terminal. This
guide is currently described for *Linux* only, but expert users should be able
to adapt it to *Windows*.

Install `git <https://github.com/git-guides/install-git>`_ and python, we
suggest to work with miniconda/anaconda. For installation see `miniconda
<https://docs.conda.io/en/latest/miniconda.html>`_ or
`anaconda <https://docs.anaconda.com/anaconda/install/linux/>`_.


Setup AvaFrame
^^^^^^^^^^^^^^

Create a new `conda environment
<https://conda.io/projects/conda/en/latest/user-guide/concepts/environments.html>`_
for AvaFrame and activate it::

  conda create --name avaframe_env
  conda activate avaframe_env

Clone the AvaFrame repository (in a directory of your choice: [YOURDIR]) and change into it::

  cd [YOURDIR]
  git clone https://github.com/avaframe/AvaFrame.git
  cd AvaFrame

AvaFrame needs following requirements as specified in the requirements.txt:

.. include:: ../requirements.txt
   :literal:

Install these with (conda-forge is needed for cmcrameri)::

  conda install -c conda-forge --file requirements.txt

or alternatively::

  pip install -r requirements.txt

Finally install avaframe::

  pip install -e .

This installs avaframe in editable mode, so every time you import avaframe the
current (local) version will be used.

Test it by starting ``python`` and do an ``import avaframe``. If no error comes
up, you are good to go.


First run
^^^^^^^^^


Follow these steps to run your first simulation:

* change into your ``AvaFrame`` directory (replace [YOURDIR]
  with your path from the previous steps)::

    cd [YOURDIR]/AvaFrame/avaframe

* run:
  ::

    python runCom1DFA.py

* a similar output should show up:
  ::

    logUtils - INFO -  Started logging at: 03.11.2020 22:42:04
    logUtils - INFO -  Also logging to: data/avaParabola/runCom1DFA.log
    runCom1DFA - INFO -  MAIN SCRIPT
    runCom1DFA - INFO -  Current avalanche: data/avaParabola
    ...


This will perform a dense flow avalanche simulation using the ``com1DFA`` module.
The results are saved to ``data/avaParabola/Outputs/com1DFA``. For a first
look at the results, got to the folder ``reports``, there you can find
a markdown report of the simulations performed including some plots of the results.

To display markdown files in a nice way use a markdown viewer of your choice.
Some other options are:

* Use the Atom editor with a markdown plugin
* If you have *pandoc* installed use this to convert it to pdf/html
* Some browsers have markdown extensions you can install easily


Workflow example
----------------

The following example should make it easier for you to find your way in AvaFrame
and setup your own AvaFrame workflow after you did the full setup.

Make sure you change to your ``AvaFrame`` direcotry by::

  cd [YOURDIR]/AvaFrame

Replace ``[YOURDIR]`` with the directory from your installation step.

Initialize project
^^^^^^^^^^^^^^^^^^

To create the folder where the input data lies and where the
output results will be saved, specify the full path to the folder
in the ``local_avaframeCfg.ini`` (which is a copy of
``avaframeCfg.ini`` that you need to create). So::

  cd avaframe
  cp avaframeCfg.ini local_avaframeCfg.ini

and edit ``local_avaframeCfg.ini`` with your favorite text editor and adjust the
variable ``avalancheDir``.

Then run ::

  python runInitializeProject.py

This will create a new directory with the input required by AvaFrame structured as described
in :ref:`moduleIn3Utils:Initialize Project`.

Input data
^^^^^^^^^^

Check the input data required by the different modules you want to use and fill
the ``Inputs/`` inside the ``[avalancheDir]`` folder from the initialize step accordingly.

For example the ``com1DFA`` module needs input as described in
:ref:`moduleCom1DFA:Input`. You can also have a look at the default setting for
the module you want to use (for example ``com1DFACfg.ini`` for module com1DFA).
If you want to use different settings, create a ``local_`` copy of the ``.ini``
file and modify the desired parameters.

More information about the configuration can be found here: :ref:`configuration:Configuration`

Building your run script
^^^^^^^^^^^^^^^^^^^^^^^^

Create your own workflow by taking the ``runOperational.py`` script as
template.

We suggest you copy it and adjust it to your liking. There are annotations in
the code that should help you to understand the structure.

A lot more examples can be found in the ``runScripts`` directory.
