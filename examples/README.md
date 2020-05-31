                      CWM 4-1-2020

This directory contains two examples using the ECS_DVRHelper class
library that implements Exterior Complex Scaling (ECS) using the
Finite Element Discrete Variable Representation numerical methods

Both examples are for H2 using the accurate potential curve fit of
Waech and Bernstein

ECS_FEMDVR_diatomic_time_dep_vibration_H2.py    
ECS_FEMDVR_diatomic_time_indep_vibration_H2.py

The Time-independent example reproduces a figure of the resonace
wave function for rotational angular momentum j = 17, which has a
centrifugal barrier to dissociation that binds metastable states.

The Time-dependent example propagates an initially Gaussian wave
packet that starts centered at a value of R just inward of the
maximum in the potential V(R) +j(j+1)/2*mu*R**2

Plotting output from both examples is written in ./Plot_Output
while .dat files are written in this directory.  Spectrum.dat
contains the eigenvalues of the ECS scaled Hamiltonian.

      =======ECS FEM-DVR class library ============

The directory DVR contains two .py files with two distinct class
libraries

DVRHelper.py
ECS_DVRHelper.py

Although the routines have the same names, the ones in ECS_DVRHelper
are specific to the implementation of ECS.  The comments in
ECS_DVRHelper.py describe the status of the implementation, which
is currently coded only for propagation on a single potential curve.

