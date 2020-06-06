.. include:: <isonum.txt>

.. role:: bolditalic
   :class: bolditalic

==========
Background
==========

These notes are an introduction to Discrete Variable Representations (DVRs) using an example that has particularly general applicability. The Finite Element Method with a Discrete Variable Representation (FEM-DVR) provides a way to solve the time-independent or time-dependent Schrödinger equation that, like all DVR methods, is more accurate and faster than finite difference. This method is one of a family of Discrete Variable Representations that are in common use today in chemistry and physics. It can be applied to problems with any potential on any interval. These notes explain the FEM-DVR method using Gauss-Lobatto quadrature, and also outline the Crank-Nicolson propagator for solving the time-dependent Schrödinger equation.

Introduction
------------

These notes describe a method for solving the Schrödinger equation for a particle moving in one dimension with coordinate :math:`x` for any potential :math:`V (x)` on any interval of :math:`x`. The variational method of course provides a way to do so, but its application generally poses a practical problem we would like to overcome: If we expand the unknown wave function in :math:`H |\Psi\rangle = E |\Psi\rangle` in a finite set of basis functions

.. math::
  |\Psi\rangle \approx \sum_{n=1}^N c_n |\varphi_n\rangle

substitute it into the Schrödinger equation, and project from the left with :math:`\langle \varphi_m |`, we come quickly to the familiar matrix representation

.. math::
  \mathbf{H}\vec{c} &= E \vec{c} \\
  \textrm{with} \quad H_{mn} &= \langle \varphi_m|\hat{T}|\varphi_n \rangle + \langle \varphi_m|\hat{V}|\varphi_n \rangle

that we can also get from the variational theorem. This is a variational basis representation of the Schrödinger equation.

To construct this matrix eigenvalue problem we need the matrix elements of both the kinetic energy and potential energy operators :math:`\hat{T}` and :math:`\hat{V}`. If we choose a basis for which the kinetic energy matrix elements are easy to evaluate, energy operators T and then try to apply it to solving this problem for various potentials, we generally find that the matrix elements are difficult to perform for many of those potential functions. The DVR is a way of getting around this problem. It is described in an extensive literature on many kinds of DVR that began in the 1980s with seminal work by John Light in the Chemistry Department at the University of Chicago and his coworkers [1]. The central idea is this: We choose a particular form of basis functions, no matter what the potential, that are constructed based on a Gaussian quadrature. Then we use that Gaussian quadrature to do every integral in the problem ­ of both the kinetic and potential energy operators. The result is that the potential energy matrix is diagonal

.. math::
  \langle \varphi_m|\hat{V}|\varphi_n\rangle = \delta_{nm}V(x_n)

where :math:`x_n` is a Gauss quadrature point. In other words, in this basis set method the potential energy matrix is always diagonal, and :bolditalic:`there are no potential energy integrals to be evaluated`. We only have to evaluate the potential on a grid of quadrature points in :math:`x`. So if we can evaluate the potential energy function as a function of the coordinates, we can reduce the Schrödinger equation directly to a matrix eigenvalue equation. To see how this works we have first to familiarize ourselves with the basics of Gaussian quadrature for performing integrals in general.

Gassian Quadrature
------------------
