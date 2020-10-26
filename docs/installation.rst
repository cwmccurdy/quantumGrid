.. highlight:: shell

============
Installation
============


Using pip
---------

To install quantumGrid using pip, run this command in your terminal:

.. code-block:: console

    $ pip install --user quantumgrid

the --user flag is so the package is installed locally and not install it on the root system, which is NOT a good idea. This is the preferred method to install packages that are in the Pypi index but the recommended way of installation is using Anaconda.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/

Using conda
-----------

Before going further into conda, you should always update your conda package system:

.. code-block:: console

    $ conda update -y conda

We would rather work with conda environments but this requires a few extra steps. First we have to pull our quantumgrid package from the Pypi index and setup a conda recipe to build the package locally. To do this we first have to get the building tools using this command:

.. code-block:: console

    $ conda install -y conda-build

Now we need to pull quantumgrid from the Pypi index and create a recipe to build our package locally for conda use:

.. code-block:: console

    $ conda skeleton pypi quantumgrid

This will create a directory called "quantumgrid" with a .yml file with our recipe to create our package.

Next we build our package:

.. code-block:: console

    $ conda-build quantumgrid

This may take a minute but once it's done, you should have a local package of quantumgrid in a bld directory (short for build) that you can then install into any conda environment you create! Here are the commands to install the quantumgrid package to a new conda environment named "qtest":

.. code-block:: console

    $ conda create --name qtest
    $ conda activate qtest
    $ conda install --use-local quantumgrid

Now you should be able to use the  quantumgrid package in your qtest conda environment! You can also run the example scripts simple by executing:

.. code-block:: console

    $ ecs_femdvr_time_indep_h2
    $ ecs_femdvr_time_dep_h2

From sources
------------

The sources for quantumGrid can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/zstreeter/quantumGrid

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/zstreeter/quantumGrid/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repo: https://github.com/zstreeter/quantumGrid
.. _tarball: https://github.com/zstreeter/quantumGrid/tarball/master
