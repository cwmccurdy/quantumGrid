.. role:: bolditalic
   :class: bolditalic

.. role:: bold
   :class: bold

.. role:: italic
   :class: italic

========
Examples
========

There are four example scripts that come with quantumGrid: :bolditalic:`ecs_femdvr_time_indep_h2` for a time independent calculation and :bolditalic:`ecs_femdvr_time_dep_h2` for a time dependent calculation and two examples that calculate vibrational states of :math:`H_2` and CO called :bolditalic:`femdvr_vib_states_h2` and :bolditalic:`femdvr_vib_states_co`, respectively. Once quantumGrid is installed in your local environment, these scripts can be called without typing "python". To see command line options for both scripts, just use the help command:

.. code-block:: console

    $ ecs_femdvr_time_dep_h2 --help

For example, to turn on plotting for the time independent example run the script with the plotting option:

.. code-block:: console

    $ ecs_femdvr_time_dep_h2 --want_to_plot=True


ecs_femdvr_time_indep_h2
------------------------

For time-independent potential, this example implements Exterior
Complex Scaling on the FEM-DVR contour.  The value of R0 and the
complex scale factor :math:`e^{I*theta}` are specified.  The representation
of the potential must be able to be evaluated on the complex part
of the contour.

Example:
   Finds all eigenvalues of complex scaled Hamiltonian and
   plots any one of them, specified by n_Plot

Args:
  1) want_to_plot (bool): Optional command that turns on plotting; default is false.

Potentials defined here:
  I) Morse potential for :math:`H_2`
  II) Bernstein fit of Kolos and Wolneiwicz potential with :math:`\frac{1}{R^6}`, :math:`\frac{1}{R^8}`, :math:`\frac{1}{R^{10}}` asymptotic behavior -- Gives near spectroscopic accuracy used in :cite:`TURNER1982127`, results there are reproduced by this code.

ecs_femdvr_time_dep_h2
------------------------

Finds all eigenvalues of complex scaled :math:`H_2` Hamiltonian
for nuclear motion plots any one of them, specified by n_Plot
Then starts a Gaussian packet in the well (e.g. with :math:`j=17`)
and follows it as it separates into a part that is bound in
the well and a part that dissociates and vanishes on the ECS
contour.

Args:
  1) number_of_time_intervals (int): First command line argument. Number of time intervals to perform the Crank-Nicolson propagation; defaults to 300.
  2) time_step (int): Second command line argument. Time step in the propagator to relax or restrict the calculation as needed; defaults to 0.1 atu.
  3) want_to_plot (bool): Optional command that turns on plotting; default is false.
  4) want_to_animate (bool): Optional command that turns on animation; default is false.

Potentials defined here:
   I) Morse potential for :math:`H_2`
   II) Bernstein fit of Kolos and Wolneiwicz potential with :math:`\frac{1}{R^6}`, :math:`\frac{1}{R^8}`, :math:`\frac{1}{R^{10}}` asymptotic behavior -- Gives near spectroscopic accuracy used in :cite:`TURNER1982127`, results there are reproduced by this code.

femdvr_vib_states_h2
--------------------

:math:`H_2` vibrational states using CI singles and doubles potential curve
from Psi4.  This potential yields a :math:`n = 0 \rightarrow 1` excitation energy
within a few wavenumbers of the value using the NIST values for
constants of diatomic molecules for :math:`H_2` in the formula
:math:`E_n = (n+\frac{1}{2})w_e - (n+\frac{1}{2})^2 w_ex_e`, which is :math:`4158 cm^{-1}`.

Shows how to
  i) Read in and interpolate a potential function known at discrete points
  ii) Use FEMDVR class to build FEM-DVR grid
  iii) Use FEMDVR class to build Hamiltonian in DVR basis
  iv) Find eigenvalues and eigenvectors of Hamiltonian
  v) Plot eigenfunctions of the Hamiltonian

Args:
  1) want_to_plot (bool): Optional command that turns on plotting; default is false.

Potential read from file used here:
   I) potcurve_CISD_H2_ccpvTZ.dat

femdvr_vib_states_co
--------------------

CO vibrational states using CI singles, doubles and triples potential curve from Psi4.
This potential gives a dissociation energy of :math:`~12.2 eV`, not very good
by comparison to the :math:`~11.1 eV` experimental value.
It yields a :math:`n = 0 \rightarrow 1` excitation energy of :math:`2207 cm^{-1}`
compared with the value using the NIST values for
constants of diatomic molecules for :math:`H_2` in the formula
:math:`E_n = (n+\frac{1}{2})w_e - (n+\frac{1}{2})^2 w_ex_e`, which is :math:`2143 cm^{-1}`
So not quite spectroscopic accuracy.

Shows how to
  i) Read in and interpolate a potential function known at discrete points
  ii) Use FEMDVR class to build FEM-DVR grid
  iii) Use FEMDVR class to build Hamiltonian in DVR basis
  iv) Find eigenvalues and eigenvectors of Hamiltonian
  v) Plot eigenfunctions of the Hamiltonian

Args:
  1) want_to_plot (bool): Optional command that turns on plotting; default is false.

Potential read from file used here:
   I) potcurve_CO_CISDT_ccpvDZ.dat

Modifying Scripts
-----------------

The actual names of these four example scripts are ECS_FEMDVR_diatomic_time_indep_vibration_H2.py, ECS_FEMDVR_diatomic_time_dep_vibration_H2.py, H2_vib_states_FEM_DVR.py, and CO_vib_states_FEM_DVR.py. If you downloaded the source package from github, then these examples are in the examples directory. If quantumgrid was installed using the conda instruction then the scripts should be in :italic:`/Path/to/Anaconda/envs/YOUR_ENVIRONMENT_NAME/lib/python3.7/site-packages/quantumgrid_examples`. If you are in a Unix environment then you can simply find them with the following command:

.. code-block:: console

    $ locate ECS_FEMDVR_diatomic_time_dep_vibration_H2.py

At any rate, once found you can modify your script however you like!

References
----------

.. bibliography:: _static/refs_examples.bib
  :style: unsrt
