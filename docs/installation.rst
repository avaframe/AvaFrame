Installation ================ This is a quick guide on how to install AvaFrame
and the required dependencies on your machine. AvaFrame is developed on **Linux
machines** (Ubuntu/Manjaro/Arch) with recent Python versions 3.8/3.9. An
(experimental) windows version is available too, but be aware this requires you
to modify some of the instructions below. These instructions assume you are
familiar with working in a terminal.

Install `git <https://github.com/git-guides/install-git>`_ and python, we 
suggest to work with miniconda/anaconda. For installation see `miniconda
<https://docs.conda.io/en/latest/miniconda.html>`_ or
`anaconda <https://docs.anaconda.com/anaconda/install/linux/>`_.


AvaFrame 
---------------

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

Install these with (conda-forge is needed for cmocean)::
 
  conda install -c conda-forge --file requirements.txt

or alternatively::

  pip install -r requirements.txt

Finally install avaframe either with conda (might require a ``conda install conda-build``)::

  conda develop .

or by using pip::

  pip install -e .

This installs avaframe in editable mode, so every time you import avaframe the
current (local) version will be used.

Test it by starting ``python`` and do an ``import avaframe``. If no error comes
up, you are good to go. 

Com1DFA Executable
------------------

The computation of the ``com1DFA`` dense flow avalanche module relies on a C++ executable.
The executable (for now 64bit linux and windows) and needed files are
available in this `git repository <https://github.com/avaframe/com1DFA_Exe>`_.
To install, change into your directory [YOURDIR] from the AvaFrame installation above and
clone the repository::
 
  cd [YOURDIR]
  git clone https://github.com/avaframe/com1DFA_Exe

Rename the exectuables according to your operating system, i.e. for Linux do::

  mv com1DFA_Exe/com1DFA_x86_64.exe com1DFA_Exe/com1DFA.exe
  mv com1DFA_Exe/SHPConv_Linux.exe com1DFA_Exe/SHPConv.exe

for Windows do::

  mv com1DFA_Exe/com1DFA_win64.exe com1DFA_Exe/com1DFA.exe
  mv com1DFA_Exe/SHPConv_win.exe com1DFA_Exe/SHPConv.exe


Go to the ``com1DFA`` directory of the AvaFrame repository from above and copy the configuration file::

  cd AvaFrame/avaframe/com1DFA
  cp com1DFACfg.ini local_com1DFACfg.ini

Open the ``local_com1DFACfg.ini`` file in you prefered text editor and change
the ``com1Exe`` variable to reflect you paths, i.e.::

  com1Exe = [YOURDIR]/com1DFA_Exe/com1DFA.exe -files [YOURDIR]/com1DFA_Exe/files/AK_Attributes


.. Attention::  We suggest to use the full path.

To test go to [YOURDIR], change into the com1DFA_Exe repository and run the
executable::

  cd [YOURDIR]
  cd com1DFA_Exe
  ./com1DFA.exe -files files/AK_Attributes/

The output should start like this::

  Setting config files directory: files/AK_Attributes/	(src/SW_Workspace.cpp:3435)
  ./com1DFA.exe -files files/AK_Attributes/ 	(src/SW_Workspace.cpp:3453)
  =================================================================
  ./com1DFA.exe
  Compiled Oct 19 2020 21:34:15
  ...

Exit by pressing ``q``
