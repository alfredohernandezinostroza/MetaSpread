Tutorial: running a simulation in interactive mode
=====

.. :contents::
    :depth: 1

This tutorial will guide you through the basic usage of MetaSpread through the interactive mode, from running a simulation to visualizing the results.

Installation
------------
If you haven't installed MetaSpread yet, use the following command:

.. code-block:: console
   
   pip install metaspread

or:

.. code-block:: console

   python -m pip install metaspread

Running metaspread
------------------

Once installed, create a folder for this tutorial project. For example, this command will create a folder called "tutorial-interactive":
.. code-block:: console

   mkdir tutorial-interactive

Then, navigate to this folder:
.. code-block:: console

   cd tutorial-interactive

Use the following command on a terminal to run MetaSpread in interactive mode:

.. code-block:: console

   python -m metaspread


Since this is the first time you run MetaSpread, it will create a file called ``simulation_configs.csv``. This file contains the default parameters for the simulation, which you can modify later.

Afterwards, the interactive interface will be launched, displaying the main menu.:

.. image:: main_menu.png

The main menu provides several options for running simulations and visualizing results. To select an option, move with the up and down arrow keys. To select an option, press enter.

First, we will focus on running a new simulation. 

Running a Simulation
--------------------

1. **Editing the parameters**: the default parameters are already set in the ``simulation_configs.csv`` file. You can edit this file to change the parameters of your simulation. To edit the file, you can use any text editor or spreadsheet software that supports CSV format. After making changes, save the file.
2. **Start a simulation**: select the ``New Simulation`` option from the main menu.
