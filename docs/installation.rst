.. highlight:: shell

============
Installation
============

Recommended
-----------

Before going further using conda, you should always update your conda package system:

.. code-block:: console

    $ conda update -y conda

First, we always recommend to create a conda environment for using our package. This is the general recommended procedure so there are no dependency issues and systemic issues that could break other parts of your computer. To create a conda environment named "DVRenv" run this command in your terminal:

.. code-block:: console

    $ conda create -n DVRenv

Now we want all the rest of the packages to only be installed into this environment so activate it before moving forward:

.. code-block:: console

    $ conda activate DVRenv

If you have Anaconda intergrated with your shell, you should see `(DVRenv)` in the front of your prompt, indicating you are now in the `DVRenv` environment. If you do not have Anaconda integrated with your shell, then run the following command and confirm you see `DVRenv` on the next line in your terminal:

.. code-block:: console

    $ echo $CONDA_DEFAULT_ENV
    $ DVRenv

Now the quantumgrid package is on the PyPI index so we need to install pip to access that index.

.. code-block:: console

   (DVRenv) $ conda install pip

(Note that this pip will only be installed in our `DVRenv` environment!)

Now we can install quantumgrid!

.. code-block:: console

   (DVRenv) $ pip install quantumgrid

Now you should be able to use the  quantumgrid package in your `DVRenv` conda environment! You can also run the example scripts simply by executing (see example directory for more details):

.. code-block:: console

    (DVRenv) $ ecs_femdvr_time_indep_h2
    (DVRenv) $ ecs_femdvr_time_indep_h2 --want_to_plot=true

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
