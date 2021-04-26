"""
    C.W. McCurdy 03/13/2020 :bold:`Example Script`
    Time-dependent Exterior Complex Scaling (ECS) FEM-DVR example
    Uses femdvr.py and potential class.

    Finite Element Method - Discrete Variable Representation (FEM-DVR)
    for 1D Schrödinger equation using Gauss-Lobatto quadrature in
    each finite element Uses class FEM_DVR(args) to construct FEM-DVR
    points, weights and Kinetic Energy

    Example:
       Finds all eigenvalues of complex scaled H2 Hamiltonian
       for nuclear motion plots any one of them, specified by n_Plot
       Then starts a Gaussian packet in the well (e.g. with :math:`j=17`)
       and follows it as it separates into a part that is bound in
       the well and a part that dissociates and vanishes on the ECS
       contour.

    Potentials defined here:
       (1) Morse potential for :math:`H_2`
       (2) Bernstein fit of Kolos and Wolneiwicz potential with :math:`\\frac{1}{R^6}`, :math:`\\frac{1}{R^8}`, :math:`\\frac{1}{R^{10}}` asymptotic behavior -- Gives near spectroscopic accuracy used in :cite:`TURNER1982127`, results there are reproduced by this code.

    References
    ----------

    .. bibliography:: _static/refs_examples.bib
      :style: unsrt
"""
# preliminaries to invoke SciPy linear algebra functions
from scipy import linalg as LA

# and NgbPy which is used to define pi, sqrt, array, .transpose etc. as
import numpy as np
import matplotlib.pyplot as plt  # import matplotlib pyplot functions
from matplotlib import animation  # for animation from same class library
import os  # functions to manipulate files and directories

import time as timeclock  # for timing parts of the calculation during debugging

# for debugging
# import sys

# sys.path.append("../")

# Importing our classes
from quantumgrid.femdvr import FEM_DVR
from quantumgrid.potential import Potential

import sys
import click


@click.option(
    "--want_to_plot_animate",
    type=click.BOOL,
    default="False",
    help="Set to True if you want to turn on plotting and animation",
)
@click.argument("time_step", type=click.FLOAT, default=0.1, required=False)
@click.argument(
    "number_of_time_intervals", type=click.INT, default=300, required=False
)
@click.command()
def main(number_of_time_intervals, time_step, want_to_plot_animate):
    #
    # ================ Make Directory for Plots if it's not there already =============
    #
    # detect the current working directory and print it
    path = os.getcwd()
    print("The current working directory is %s" % path)
    # define the name of the directory to be created
    Plot_Output = path + "/Plot_Output"

    if want_to_plot_animate is True:
        if os.path.exists(Plot_Output):
            print(
                "Directory for wave function plots already exists", Plot_Output
            )
        else:
            print(
                "Attempting to create directory for wave function plots ",
                Plot_Output,
            )
            try:
                os.mkdir(Plot_Output)
            except OSError:
                print("Creation of the directory %s failed" % Plot_Output)
            else:
                print("Successfully created the directory %s " % Plot_Output)
    # =====================================FEM_DVR===================================
    #  Set up the FEM DVR grid given only the Finite Element boundaries and order
    #  of Gauss Lobatto quadrature,  and compute the Kinetic Energy matrix for
    #  this FEM DVR for Mass = mass set in call (atomic units).
    #
    # Here is where the reduced mass is set
    # H_Mass = 1.007825032 #H atom atomic mass
    # H_Mass = 1.00727647  #proton atomic mass
    # standard atomic weight natural abundance, seems to be what Turner, McCurdy used
    H_Mass = 1.0078
    Daltons_to_eMass = 1822.89
    mu = (H_Mass / 2.0) * Daltons_to_eMass
    print("reduced mass mu = ", mu)
    bohr_to_Angstrom = 0.529177
    Hartree_to_eV = 27.211386245988  # NIST ref
    eV_to_wavennumber = 8065.54393734921  # NIST ref on constants + conversions
    # value from NIST ref on constants + conversions
    Hartree_to_wavenumber = 2.1947463136320e5
    HartreeToKelvin = 315773
    atu_to_fs = 24.18884326509 / 1000
    #  Set up the FEM-DVR grid
    n_order = 35
    FEM_boundaries = [
        0.4,
        1.0,
        2.0,
        3.0,
        4.0,
        5.0,
        6.0,
        7.0,
        8.0,
        10.0,
        11.0,
        12.0,
        15.0,
        20.0,
        25.0,
    ]
    # parameters for Fig 2 in Turner, McCurdy paper
    # Julia Turner and C. William McCurdy, Chemical Physics 71(1982) 127-133
    scale_factor = np.exp(1j * 20.0 * np.pi / 180.0)
    R0 = 10.0  # potential must be analytic for R >= R0
    fem_dvr = FEM_DVR(
        n_order,
        FEM_boundaries,
        Mass=mu,
        Complex_scale=scale_factor,
        R0_scale=R0,
    )
    print("\nFEM-DVR basis of ", fem_dvr.nbas, " functions")

    pertubation = Potential()
    # ==================================================================================
    #  Plot potential on the DVR grid points on which the wavefunction is defined
    print("\n Plot potential ")
    time = 0.0
    x_Plot = []
    pot_Plot = []
    for j in range(0, fem_dvr.nbas):
        x_Plot.append(np.real(fem_dvr.x_pts[j + 1]))
        pot_Plot.append(
            np.real(pertubation.V_Bernstein(fem_dvr.x_pts[j + 1], time))
        )

    if want_to_plot_animate is True:
        plt.suptitle(
            "V(x) at DVR basis function nodes", fontsize=14, fontweight="bold"
        )
        string = "V"
        plt.plot(x_Plot, pot_Plot, "ro", label=string)
        plt.plot(x_Plot, pot_Plot, "-b")
        plt.legend(loc="best")
        plt.xlabel(" x ", fontsize=14)
        plt.ylabel("V", fontsize=14)
        print(
            "\n Running from terminal, close figure window to proceed and make .pdf file of figure"
        )
        #   Insert limits if necessary
        #   Generally comment this logic.  Here I am matching the Turner McCurdy Figure 2
        # CWM: need to use float() to get plt.xlim to work to set x limits
        ymax = float(0.05)
        plt.ylim([-0.18, ymax])
        # save plot to .pdf file
        plt.savefig(
            "Plot_Output/" + "Plot_potential" + ".pdf", transparent=False
        )
        plt.show()
    #
    # =============Build Hamiltonian (at t=0 if time-dependent)=================================
    #     Pass name of potential function explicitly here
    time = 0.0
    H_mat = fem_dvr.Hamiltonian(pertubation.vectorized_V_Bernstein, time)
    print("\n Completed construction of Hamiltonian at t = 0")
    # ====================================================================================
    #
    # Find all the eigenvalues of the Hamiltonian so we can compare with known bound state energies
    # or make a plot of the spectrum -- For a time-independent Hamiltonian example here
    #
    # EigenVals = LA.eigvals(H_mat) # eigenvalues only for general matrix.  LA.eigvalsh for Hermitian
    print(
        "Calculating ",
        fem_dvr.nbas,
        " eigenvalues and eigenvectors for plotting eigenfunctions",
    )
    EigenVals, EigenVecs = LA.eig(H_mat, right=True, homogeneous_eigvals=True)
    print("after LA.eig()")
    #
    n_energy = fem_dvr.nbas
    file_opened = open("Spectrum_ECS.dat", "w")
    print("EigenVals shape ", EigenVals.shape)
    for i in range(0, n_energy):
        print("E( ", i, ") =   ", EigenVals[0, i], " hartrees")
        print(
            np.real(EigenVals[0, i]),
            "  ",
            np.imag(EigenVals[0, i]),
            file=file_opened,
        )
    # ====================================================================================
    #
    # Extract the n_Plot'th eigenfunction for plotting
    #
    # pick one of the bound states of Morse Potential to plot
    # numbering can depend on numpy and python installation that determines
    # behavior of the linear algebra routines.
    n_Plot = (
        n_energy - 1
    )  # This is generally the highest energy continuum eigenvalue
    n_Plot = 426
    wfcnPlot = []
    for j in range(0, fem_dvr.nbas):
        wfcnPlot.append(EigenVecs[j, n_Plot])
    #
    # Normalize  wave function from diagonalization
    # using integration of square on the contour
    # following the original Rescigno McCurdy idea for partial widths
    # from the paper
    # "Normalization of resonance wave functions and the calculation of resonance widths"
    #  Rescigno, McCurdy, Phys Rev A 34,1882 (1986)
    norm_squared = 0.0
    for j in range(0, fem_dvr.nbas):
        norm_squared = norm_squared + (wfcnPlot[j]) ** 2
    wfcnPlot = wfcnPlot / np.sqrt(norm_squared)
    norm_squared = 0.0
    gamma_residue = 0.0
    # background momentum defined with Re(Eres)
    k_momentum = np.sqrt(2.0 * mu * np.real(EigenVals[0, n_Plot]))
    # k_momentum = np.real(np.sqrt(2.0*mu*EigenVals[0,n_Plot]))
    for j in range(0, fem_dvr.nbas):
        norm_squared = norm_squared + (wfcnPlot[j]) ** 2
        if fem_dvr.x_pts[j + 1] < 5.8:
            free_wave = (2.0 * np.sqrt(mu / k_momentum)) * np.sin(
                k_momentum * fem_dvr.x_pts[j + 1]
            )
            gamma_residue = gamma_residue + wfcnPlot[
                j
            ] * pertubation.V_Bernstein(
                fem_dvr.x_pts[j + 1], time
            ) * free_wave * np.sqrt(
                fem_dvr.w_pts[j + 1]
            )
    print("Complex symmetric inner product (psi|psi) is being used")

    if want_to_plot_animate is True:
        print(
            "Norm of wave function from int psi^2 on contour being plotted is ",
            np.sqrt(norm_squared),
        )

    print(" For this state the asymptotic value of k = ", k_momentum)
    print(
        "gamma from int = ",
        gamma_residue,
        " |gamma|^2 = ",
        np.abs(gamma_residue) ** 2,
    )
    # Plot wave function -- It must be type np.complex
    Cinitial = np.zeros((fem_dvr.nbas), dtype=np.complex)
    wfcnInitialPlot = np.zeros((fem_dvr.nbas), dtype=np.complex)
    for j in range(0, fem_dvr.nbas):
        Cinitial[j] = wfcnPlot[j]

    if want_to_plot_animate is True:
        #
        # plot n_Plot'th eigenfunction
        #
        print(
            "\n Plot Hamiltonian eigenfunction number ",
            n_Plot,
            " with energy ",
            EigenVals[0, n_Plot],
        )
        tau = atu_to_fs / (-2.0 * np.imag(EigenVals[0, n_Plot]))
        print("  Lifetime tau = 1/Gamma = ", tau, " fs")
        number_string = str(n_Plot)
        title = "Wavefunction number = " + number_string
        wfn_plot_points = 2000
        x_Plot_array, Psi_plot_array = fem_dvr.Plot_Psi(
            Cinitial,
            plot_title_string=title,
            N_plot_points=wfn_plot_points,
            make_plot=True,
        )
        #
        #  Make data file for n_Plot'th eigenfunction
        #  Psi and the integrand of the residue factors, gamma, of the free-free matrix element
        #  of the full Green's function, as per
        #
        filename = "wavefunction" + number_string + ".dat"
        file_opened = open(filename, "w")
        print_points = len(x_Plot_array)
        print("x_Plot_array shape ", print_points)
        for i in range(print_points):
            free_wave = (2.0 * np.sqrt(mu / k_momentum)) * np.sin(
                k_momentum * x_Plot_array[i]
            )
            # for partial width gamma
            integrand = (
                Psi_plot_array[i]
                * pertubation.V_Bernstein(x_Plot_array[i], time)
                * free_wave
            )
            print(
                np.real(x_Plot_array[i]),
                "  ",
                np.imag(x_Plot_array[i]),
                "  ",
                np.real(Psi_plot_array[i]),
                "  ",
                np.imag(Psi_plot_array[i]),
                "  ",
                np.real(integrand),
                np.imag(integrand),
                file=file_opened,
            )
    #
    # ================# Initialize wave function at t = 0 ================================
    #   It must be type np.complex
    Cinitial = np.zeros((fem_dvr.nbas), dtype=np.complex)
    wfcnInitialPlot = np.zeros((fem_dvr.nbas), dtype=np.complex)
    for j in range(0, fem_dvr.nbas):
        #
        #  Displaced Gaussian packet in well. Displace to larger R because displacing
        #  to R < Req produces significant amplitude in the continuum
        #
        alpha = 40.0
        shift = 4.0
        Cinitial[j] = (
            np.sqrt(np.sqrt(2.0 * alpha / np.pi))
            * np.exp(-alpha * (fem_dvr.x_pts[j + 1] - shift) ** 2)
            * np.sqrt(fem_dvr.w_pts[j + 1])
        )
        wfcnInitialPlot[j] = Cinitial[j]
    #
    # =====================Propagation from t = tinitial to tfinal ============================
    # initialize arrays for storing the information for all the frames in the animation of
    # the propagation of the wave function
    times_array = []
    x_Plot_time_array = []
    Psi_plot_time_array = []
    tinitial = 0.0
    # specify final time and specify number of intervals for plotting
    #  commented info for Morse oscillator revivals, in atu
    # a = 1.0277
    # d = 0.1746
    # omega = np.sqrt((2.e0/mu)*d*a**2)
    # T_revival = 4*np.pi*mu/a**2 # for the morse oscillator
    # number_of_time_intervals = 25*np.int((tfinal-tinitial)*omega/(2.0*np.pi))
    # print("T_revival = ",T_revival," atomic time units ",T_revival*24.1888/1000.0," femtoseconds")
    tfinal = 5000  # specified in atomic time units
    print(
        "\nOverall propagation will be from ",
        tinitial,
        " to ",
        tfinal,
        " atu in ",
        number_of_time_intervals,
        " intervals",
    )

    if want_to_plot_animate is True:
        #
        # plot initial wave packet
        #
        print("\n Plot initial wave function at t = ", tinitial)
        number_string = str(0.0)
        title = "Wavefunction at t = " + number_string
        x_Plot_array, Psi_plot_array = fem_dvr.Plot_Psi(
            Cinitial,
            plot_title_string=title,
            N_plot_points=1250,
            make_plot=True,
        )
        times_array.append(tinitial)
        x_Plot_time_array.append(x_Plot_array)
        Psi_plot_time_array.append(Psi_plot_array)
    #
    #  Loop over number_of_time_intervals intervals that make up t = 0 to tfinal
    #
    t_interval = (tfinal - tinitial) / number_of_time_intervals
    # for H2 reduced mass, Deltat = 0.05 in Morse oscillator seems to work for 10s of fs
    # for D2 reduced mass, Deltat = 0.1 in Morse oscillator seems to work for at least 2 ps
    # based on comparison with time steps 5 to 10 times smaller
    time_step = 0.1  # atomic time units
    # time_step = 0.05  # atomic time units
    # time_step = 0.015  # atomic time units
    for i_time_interval in range(0, number_of_time_intervals):
        t_start = tinitial + i_time_interval * t_interval
        t_finish = tinitial + (i_time_interval + 1) * t_interval
        N_time_steps = np.int(t_interval / time_step)
        Delta_t = t_interval / np.float(N_time_steps)
        #
        # Check norm of initial wavefunction on real part of ECS Contour
        #
        # find  FEM DVR function at the last real DVR point
        first_ECS_element = fem_dvr.i_elem_scale
        last_real_dvr_fcn = 1
        for i in range(0, fem_dvr.nbas):
            if fem_dvr.x_pts[i + 1] < FEM_boundaries[first_ECS_element]:
                last_real_dvr_fcn = last_real_dvr_fcn + 1
        # the point last_real_dvr_fcn still has complex bridging function weight.
        # small error for integral (0,R0) may result.
        norm_initial = 0
        for j in range(0, last_real_dvr_fcn):
            norm_initial = norm_initial + np.abs(Cinitial[j]) ** 2
        print(
            "\nNorm of starting  wave function on real part of ECS contour",
            norm_initial,
        )
        #
        #  use the propagator from DVR()
        #
        clock_start = timeclock.time()
        Ctfinal = fem_dvr.Crank_Nicolson(
            t_start,
            t_finish,
            N_time_steps,
            Cinitial,
            pertubation.vectorized_V_Bernstein,
            Is_H_time_dependent=False,
        )
        clock_finish = timeclock.time()
        print(
            "Time for ",
            N_time_steps,
            " Crank-Nicolson steps was ",
            clock_finish - clock_start,
            " secs.",
        )
        #
        # Check norms of initial and final wavefunctions
        #
        norm_final = 0
        for j in range(0, last_real_dvr_fcn):
            # should be change to real part of contour.
            norm_final = norm_final + np.abs(Ctfinal[j]) ** 2
        print(
            "Norm of final wave function on real part of ECS contour",
            norm_final,
        )

        if want_to_plot_animate is True:
            #
            # Plot of packet at end of each interval using Plot_Psi from DVR()
            #
            print(
                "Plot function propagated to t = ",
                t_finish,
                " atomic time units ",
                t_finish * atu_to_fs,
                " fs",
            )
            number_string = str(t_finish)
            title = "Wavefunction at t = " + number_string
            x_Plot_array, Psi_plot_array = fem_dvr.Plot_Psi(
                Ctfinal,
                plot_title_string=title,
                N_plot_points=1250,
                make_plot=False,
            )
            times_array.append(t_finish)
            x_Plot_time_array.append(x_Plot_array)
            Psi_plot_time_array.append(Psi_plot_array)
        #
        #   reset initial packet as final packet for this interval
        for i in range(0, fem_dvr.nbas):
            Cinitial[i] = Ctfinal[i]
    #
    if want_to_plot_animate is True:
        x_Plot_array, Psi_plot_array = fem_dvr.Plot_Psi(
            Ctfinal, plot_title_string=title, N_plot_points=1250, make_plot=True
        )
        print(
            "Final frame at tfinal = ",
            tfinal,
            "atomic time units  is showing -- ready to prepare animation ",
        )
        print("\n **** close the plot window to proceed ****")
    # ==============================================================================
    # # initialization function: plot the background of each frame
    # ==============================================================================

    def init():
        line.set_data([], [])
        time_text.set_text("")
        ax.set_xlabel(" x (bohr) ", fontsize=16, fontweight="bold")
        ax.set_ylabel(" Psi(t): Re, Im & Abs ", fontsize=16, fontweight="bold")
        fig.suptitle(
            "Wave Packet in H2 Potential", fontsize=16, fontweight="bold"
        )
        # put in a line at the value Phi = 0
        ax.plot([x_Plot[0], x_Plot[len(x_Plot) - 1]], [0, 0], "k")
        return line, time_text

    nframes = number_of_time_intervals

    def animate(i):
        # ==============================================================================
        #
        #  commands to specify changing elements of the animation
        #  We precomputed all the lines at all the times so that this
        #  wouldn't recompute everything at each frame, which is very slow...
        #
        # ==============================================================================
        time_string = str(
            times_array[i] * 24.189 / 1000.0
        )  # string for a plot label
        re_array = np.real(Psi_plot_time_array[i])
        im_array = np.imag(Psi_plot_time_array[i])
        abs_array = np.abs(Psi_plot_time_array[i])
        line1.set_data(x_Plot_array, re_array)
        line2.set_data(x_Plot_array, im_array)
        line3.set_data(x_Plot_array, abs_array)
        time_text.set_text("time = " + time_string + " fs")
        return (line1, line2, line3, time_text)

    if want_to_plot_animate is True:
        # ==============================================================================
        #  Now plot the animation
        # ==============================================================================
        # reinitialize the figure for the next plot which is the animation
        fig = plt.figure()
        ymax = 1.25
        xmin = x_Plot[0]
        xmax = 20.0
        ax = fig.add_subplot(
            111, autoscale_on=False, xlim=(xmin, xmax), ylim=(-ymax, +ymax)
        )
        (line,) = ax.plot([], [], "-r", lw=2)
        (line1,) = ax.plot([], [], "-r", lw=2)
        (line2,) = ax.plot([], [], "-b", lw=2)
        (line3,) = ax.plot([], [], "k", lw=2)
        # define the object that will be the time printed on each frame
        time_text = ax.text(0.02, 0.95, "", transform=ax.transAxes)

        # ==============================================================================
        # call the animator.  blit=True means only re-draw the parts that have changed.
        # ==============================================================================
        anim = animation.FuncAnimation(
            fig,
            animate,
            init_func=init,
            frames=nframes + 1,
            interval=200,
            blit=True,
            repeat=True,
        )
        # ==============================================================================
        #  show the animation
        # ==============================================================================
        anim.save("Plot_Output/H2_wavepacket.mp4")
        plt.show()
        print("done")

    if want_to_plot_animate is False:
        print(
            "\n\nSet the command line option want_to_plot_animate=True to see figures,"
            " animation, and create plotting directory.\n\n"
        )


if __name__ == "__main__":
    main()
