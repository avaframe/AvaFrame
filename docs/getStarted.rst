Quick start
================


Follow these steps to run your first simulation:

* go to :ref:`installation:Installation` and follow the instructions to install
  *AvaFrame* and the necessary *com1DFA* executable.
* change into your ``AvaFrame`` directory (replace [YOURDIR]
  with your path from the first step)::

    cd [YOURDIR]/AvaFrame/avaframe

* run:
  ::

    python runCom1DFA.py

* a similar output should show up:
  ::

    logUtils - INFO -  Started logging at: 03.11.2020 22:42:04
    logUtils - INFO -  Also logging to: data/avaHelix/runCom1DFA.log
    runCom1DFA - INFO -  MAIN SCRIPT
    runCom1DFA - INFO -  Current avalanche: data/avaHelix
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
