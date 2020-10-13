This directory contains two examples using the FEM_DVR class
library that implements Exterior Complex Scaling (ECS) using the
Finite Element Discrete Variable Representation numerical methods and
also using the Potential class.

Both examples for distribution are for H2 using the accurate potential
curve fit of Waech and Bernstein (referenced in the Turner-McCurdy
paper below) T. G. Waech R.B. Bernstein, J. Chem. Phys. 46 (1967)
4905.

ECS_FEMDVR_diatomic_time_dep_vibration_H2.py
ECS_FEMDVR_diatomic_time_indep_vibration_H2.py

The Time-independent example reproduces a Figure 2 of Julia Turner
and C. William McCurdy, Chemical Physics 71(1982) 127-133 of the
resonace wave function for rotational angular momentum j = 17, which
has a centrifugal barrier to dissociation that binds metastable
states.  The complex resonance energy is  E_res = (0.004044878419994
-0.000219496448j)  hartrees.

The Time-dependent example propagates an initially Gaussian wave
packet that starts centered at a value of R just inward of the
maximum in the potential V(R) +j(j+1)/2*\mu*R^2

Plotting output from both examples is written in ./Plot_Output
while .dat files are written in this directory.  Spectrum.dat
contains the eigenvalues of the ECS scaled Hamiltonian.
