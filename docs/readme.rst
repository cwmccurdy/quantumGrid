About
=====
quantumGrid is a package for solving a 1-D Schrödinger equation for an
arbitrary potential on any interval. The heart of this package is using
a Finite Element Method with a Discrete Variable Representation
(FEM-DVR) grid to solve the time-dependent or time-independent
Schrödinger equation. This grid provides a compact supported foundation
for numerically accurate integration and also allows for a natural
application of outgoing scattering boundary conditions by adding a complex
tail as the last finite element of the FEM-DVR grid, called *exterior
complex scaling* (ECS). Therefore, this grid can be applied to
scattering problems where the resonances become square integrable
under this complex rotation of the Hamiltonian.


Motivation
==========

This python package was created for a graduate course in time-dependent
quantum mechanics at UC Davis. Given the generality and usefulness of a
Finite Element Method - Discrete Variable Representation (FEM-DVR) grid
for solving the Schrödinger equation and simple scattering problems, we
wanted to go open source and provide this numerical tool for others in
the name of science!

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
   <img src="_static/images/Bills_pic.jpg" width="100px;" alt=""/>

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

   </tr>

.. raw:: html

   </table>

This project follows the
`all-contributors <https://github.com/all-contributors/all-contributors>`__
specification. Contributions of any kind welcome!
