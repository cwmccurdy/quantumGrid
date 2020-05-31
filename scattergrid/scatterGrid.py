"""
Finite Element Method -- Discrete Variable Representation (DVR) for
  1D Schroedinger equation using Gauss-Lobatto quadrature
     C. William McCurdy and Giuseppe Barbalinardo -- UC Davis

  Includes time propagation with Crank-Nicolson propagator
  and routines for dynamics on two coupled potential curves.

  November 2019: time-dependent potentials allowed for two-state propagation
                 vectorized logic in Hamiltonian and Potential builds

  March 2020:  Exterior Complex Scaling implemented for single state portion
               both time-independent and time-dependent calculations
               Vectorized logic for diagonals of Hamiltonian can be commented and
               alternate logic uncommented for potential functions that don't
               vectorize correctly.  2-state routines are not yet implemented for ECS

"""
