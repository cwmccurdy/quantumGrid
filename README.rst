.. raw:: html

   <h1 align="center">
     💥quantumGrid💥
   </h1>

.. raw:: html

   <p align="center">
     <img src=docs/_static/images/scatteredWave.png width="500px;" height="400px;" alt=""/>
   </p>

brought to you by the AMO theory group at Berkeley National Lab

Note:
  The image above is for a scattered wave on a 2D FEM-DVR grid and created using POV-Ray. A 2D example is included with this package but uses Python for plotting


.. image:: https://img.shields.io/pypi/v/quantumgrid.svg
        :target: https://pypi.python.org/pypi/quantumgrid

.. image:: https://travis-ci.com/zstreeter/quantumGrid.svg?branch=master
        :target: https://travis-ci.com/zstreeter/quantumGrid

.. image:: https://readthedocs.org/projects/quantumgrid/badge/?version=latest
        :target: https://quantumgrid.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/zstreeter/quantumGrid/shield.svg
     :target: https://pyup.io/repos/github/zstreeter/quantumGrid/
     :alt: Updates

.. image:: https://pyup.io/repos/github/zstreeter/quantumGrid/python-3-shield.svg
     :target: https://pyup.io/repos/github/zstreeter/quantumGrid/
     :alt: Python 3

Table of Contents
=================

-  `About <#about>`__
-  `Motivation <#motivation>`__
-  `Installation <#installation>`__
-  `Features <#features>`__
-  `Contributors <#contributors>`__
-  `License <#license>`__

About
=====

Exterior Complex Scaled Finite-Element Element Discrete Variable
Representation grid for general physics problems. In other words,
quantumGrid is a package for solving a 1-D Schrödinger equation
for an arbitrary potential.

Motivation
==========

This python package was created for a graduate course in time-dependent
quantum mechanics at UC Davis. Given the ease of programming in python,
generality and usefulness of a Finite Element Method - Discrete
Variable Representation (FEM-DVR) grid for solving the Schrödinger
equation and simple scattering problems, we wanted to go open source
and provide this numerical tool for others in the name of science!

Features
========

.. raw:: html

   <p align="center">
     <img src=docs/_static/images/DVR.png width="800px;" height="400px;" alt=""/>
   </p>

Following the documentation link below, a few toy problems can be found
in the examples directory. There you will find example scripts
implementing time-independent and time-dependent calculations of molecular
hydrogen and two scripts for finding the vibrational states of molecular
hydroden and carbon dioxide.

Documentation
==============
For more details in using quantumGrid, checkout our manual here:
https://quantumgrid.readthedocs.io.

Contributors ✨
===============

Thanks goes to these wonderful people (`emoji
key <https://allcontributors.org/docs/en/emoji-key>`__):

.. raw:: html

   <table>

.. raw:: html

   <tr>

.. raw:: html

   </tr>

.. raw:: html

   <td align="center">
   <a href="https://chemistry.ucdavis.edu/people/william-mccurdy">
   <img src="docs/_static/images/Bills_pic.jpg" width="100px;" alt=""/>

Willaim (Bill) McCurdy 💻 🚧 📖

.. raw:: html

   <td align="center">
   <a href="https://www.linkedin.com/in/zachary-streeter-44a323102/">
   <img src="https://avatars0.githubusercontent.com/u/15461329?v=4" width="100px;" alt=""/>

Zachary Streeter💻 🚧 📖

.. raw:: html

   </td>

.. raw:: html

   <td align="center">
   <a href="http://giuseppe.barbalinardo.com">
   <img src="https://avatars2.githubusercontent.com/u/6192485?v=4" width="100px;" alt=""/>

Giuseppe Barbalinardo💻

.. raw:: html

   </td>

.. raw:: html

   </table>

This project follows the
`all-contributors <https://github.com/all-contributors/all-contributors>`__
specification. Contributions of any kind welcome!


Credits
-------

* Free software: MIT license

* This package template was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
