.. role:: bolditalic
   :class: bolditalic

.. role:: bold
   :class: bold

.. role:: italic
   :class: italic

========
Examples
========

There are two example scripts that come with quantumGrid: :bolditalic:`ecs_femdvr_time_indep_h2` for a time independent calculation and :bolditalic:`ecs_femdvr_time_dep_h2` for a time dependent calculation. Once quantumGrid is installed in your local environment, these scripts can be called without typing "python".

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

Modifying Scripts
-----------------

The actual names of these two example scripts are ECS_FEMDVR_diatomic_time_indep_vibration_H2.py and ECS_FEMDVR_diatomic_time_dep_vibration_H2.py. If quantumgrid was installed using the conda instruction then the scripts should be in :italic:`/Anaconda/envs/YOUR_ENVIRONMENT_NAME/lib/python3.7/site-packages/quantumgrid`. If you are in a Unix environment then you can simply find them with the following command:

.. code-block:: console

    $ locate ECS_FEMDVR_diatomic_time_dep_vibration_H2

Once found you can modify your script however you like!

References
----------

.. bibliography:: _static/refs_examples.bib
  :style: unsrt
