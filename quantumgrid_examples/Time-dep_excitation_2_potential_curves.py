#%%
"""                  Chem 210A/B  C.W. McCurdy

     Finite Element Method - Discrete Variable Representation (FEM-DVR)
     for 1D Schrödinger equation using Gauss-Lobatto quadrature in
     each finite element.  Uses class FEM_DVR() to construct FEM-DVR
     points, weights and Kinetic Energy and propagate wave packets.

     In this example:
        (1) Excitation from one potential curve (electronic state) to another in a diatomic molecule by a finite pulse.
        (2) The dipole matrix element between the two states is assumed to be constant as a function of internuclear distance.
        (3) Potential curves are Morse oscillator functions with different well depths and shapes shifted by  0.15 hartrees
        (4) Reduced mass is reduced mass of H2.  Note that a denser DVR grid may be necessary for heavier masses.
        (5) Pulse has Sin^2 envelope with 3 femtosecond duration and is centered at 0.2 hartrees

     Numerical parameters of DVR and propagation that control convergence of the
      numerical methods are specified by the variables
         * n_order
         * FEM_boundaries
         * time_step


     Output includes an animation of the wave packet, plots of the potentials and initial
      wave packet, and text output of grid parameters, Hamiltonian eigenvalues, and properties
      of wave packets during the propagation at the plotting intervals in time.

     Plot output, including an mp4 file of the animation is placed in the directory Plot_Output/

     Uses the basic FEM-DVR functions in the FEM_DVR() class and also the functions
       specific to the case of nuclear motion on two potential surfaces
           Hamiltonian_Two_States    <---- constructs the two Hamiltonian matrices and coupling
           Crank_Nicolson_Two_States <----- propagates the two component wave function with coupling between them
           Plot_Psi_Two_States       <---- plots both components of the wave function (on the two potentials)

"""
#%%
# preliminaries to invoke SciPy linear algebra functions
from scipy import linalg as LA

# and NumPy which is used to define pi, sqrt, array, .transpose etc. as
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
    massH = 1836.0
    mu = massH / 2.0
    n_order = 30  # 40th order with finite elements of 1 or 2 bohr OK for light masses here
    FEM_boundaries = [0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0]
    #
    # constants and conversions
    #
    Daltons_to_eMass = 1822.89
    bohr_to_Angstrom = 0.529177
    Hartree_to_wavenumber = 8065.54 * 27.211
    Hartree_to_eV = 27.211
    atu_to_fs = 24.18884 / 1000.0
    c_light = 137.035999084  # speed of light in atomic units
    a0_to_m = bohr_to_Angstrom * 10 ** (-10)
    I_0 = (
        3.5095 * 10 ** 16
    )  # conversion for intensity in W/cm^2 I =I_0*(A_0*omega/c_light)**2
    #
    print("Mass = ", mu)
    fem_dvr = FEM_DVR(n_order, FEM_boundaries, Mass=mu)
    print("\nFEM-DVR basis of ", fem_dvr.nbas, " functions")
    #
    #  pulse parameters for the dipole coupling between the two
    #  electronic states
    #
    T_pulse = 3.0  # set in femtoseconds
    omega = 0.2  # central frequency 0.113903 hartrees  is 400 nm
    A_0 = 10.0  # strength of A field 6.422 atomic units is 10^12 W/cm^2 for 400 nm
    print("\n\n*** Radiation Pulse ***")
    print(" Duration = ", T_pulse, " fs,  ")
    print(
        " Intensity of pulse  = ",
        "%15.3e" % (I_0 * (A_0 * omega / c_light) ** 2),
        " W/cm^2",
    )
    print(
        " Central energy = ",
        omega * Hartree_to_eV,
        " eV,  lambda = ",
        (10 ** 9) * a0_to_m * 2 * np.pi * c_light / omega,
        " nm",
    )

    def E_field(time):
        #
        #  Definition of E in - mu E(t) radiative coupling between electronic states
        #  One pulse here, followed by A = 0, durations and other parameters set in main program
        #
        if 0 <= time * atu_to_fs <= T_pulse:
            T = T_pulse / atu_to_fs  # pulse duration  converted from fs
            field = (
                (A_0 * omega / c_light)
                * (np.sin(np.pi * time / T) ** 2)
                * np.cos(omega * time - omega * T / 2)
            )
        else:
            field = 0.0
        return field

    def Vcoupling(x, time):
        #
        #  mu_dipole * E(t)
        #
        field = E_field(time)
        mu_dipole = 1.0  # assumed value of dipole matrix element
        coupling = mu_dipole * field
        return coupling

    def Vzero(x, time):
        Vcoupzero = 0.0
        return Vcoupzero

    pertubation = Potential()

    # =============Build Hamiltonian without coupling to pick an initial state=========
    #     Pass name of potential function explicitly here
    time = 0.0
    H_mat = fem_dvr.Hamiltonian_Two_States(
        pertubation.V_morse_1, pertubation.V_morse_2, Vzero, time
    )
    print("\n Completed construction of Hamiltonian at t = 0")
    # ====================================================================================
    #
    # Find all the eigenvalues of the Hamiltonian so we can compare with known bound state energies
    # or make a plot of the spectrum -- For a time-independent Hamiltonian example here
    #
    EigenVals_diabatic = LA.eigvalsh(H_mat)
    # ====================================================================================
    #
    # Print the lowest few eigenvalues as a check, and help pick an initial state in one diabatic well
    print(
        "\nEigenvalues of the uncoupled Hamiltonians for the two potential curves"
    )
    print("n      E_uncoupled")
    for i in range(0, 50):
        print(i, "    ", EigenVals_diabatic[i])
    #
    # ===========================================================================================
    #  Plot potential on the DVR grid points on which the wavefunction is defined
    time = (1 / 2) * T_pulse / atu_to_fs  # Maximum of pulse
    #
    print("\n Plot potential ")
    x_Plot = []
    pot1_Plot = []
    pot2_Plot = []
    coup_Plot = []
    for j in range(0, fem_dvr.nbas):
        x_Plot.append(fem_dvr.x_pts[j + 1])
        pot1_Plot.append(pertubation.V_morse_1(fem_dvr.x_pts[j + 1], time))
        pot2_Plot.append(pertubation.V_morse_2(fem_dvr.x_pts[j + 1], time))
        coup_Plot.append(10.0 * Vcoupling(fem_dvr.x_pts[j + 1], time))
    if want_to_plot_animate is True:
        plt.suptitle(
            "Uncoupled potentials and max value of coupling",
            fontsize=14,
            fontweight="bold",
        )
        string = "V1"
        # plt.plot(x_Plot,pot1_Plot,'ro',label=string) # plot points only
        plt.plot(x_Plot, pot1_Plot, "-r", label=string)
        string = "V2"
        plt.plot(x_Plot, pot2_Plot, "-b", label=string)
        string = "Max Vcoup x 10"
        plt.plot(x_Plot, coup_Plot, "-g", label=string)
        plt.legend(loc="best")
        plt.xlabel(" x (bohr)", fontsize=14)
        plt.ylabel("V (hartrees)", fontsize=14)
        print(
            "\n Running from terminal, close figure window to proceed and make .pdf file of figure"
        )
        #   Insert limits if necessary
        # xmax = float(rmax)  # CWM: need to use float() to get plt.xlim to work to set x limits
        # plt.xlim([0,xmax])
        plt.ylim([-0.2, 0.2])
        # number_string = str(a)
        plt.savefig(
            "Plot_Output/" + "Plot_two_state_potentials" + ".pdf",
            transparent=False,
        )
        plt.show()
    #
    # ====================================================================================
    #
    # Extract the n_Plot'th eigenfunction of uncoupled Hamiltonians for plotting
    #
    #
    np1 = 20
    number_of_eigenvectors = np1
    #  Initial wave function chosen here -- change n_Plot to change it
    n_Plot = 0  # pick one of the bound states of Potential to plot and use as initial state below
    print("Calculating ", np1, " eigenvectors for plotting eigenfunctions")
    EigenVals = []
    EigenVals, EigenVecs = LA.eigh(H_mat, eigvals=(0, number_of_eigenvectors))
    wfcnPlot = []
    for j in range(0, 2 * fem_dvr.nbas):
        wfcnPlot.append(EigenVecs[j, n_Plot])
    #
    # normalize 2 component wave function from diagonalization
    # so that it has the right overall normalization
    #
    norm_squared = 0.0
    for j in range(0, 2 * fem_dvr.nbas):
        norm_squared = norm_squared + np.abs(wfcnPlot[j]) ** 2
    wfcnPlot = wfcnPlot / np.sqrt(norm_squared)
    norm_squared = 0.0
    for j in range(0, 2 * fem_dvr.nbas):
        norm_squared = norm_squared + np.abs(wfcnPlot[j]) ** 2
    print("Norm of wave function being plotted is ", np.sqrt(norm_squared))
    # now extract the two components
    wfcnPlot1 = []
    wfcnPlot2 = []
    for j in range(0, fem_dvr.nbas):
        wfcnPlot1.append(
            wfcnPlot[j,]
        )
        wfcnPlot2.append(wfcnPlot[j + fem_dvr.nbas])
    #
    # Uncomment to both components of the eigenfunction.  For uncoupled Hamiltonians
    # only one component will be nonzero for each eigenfunction.
    #
    # x_Plot_array, Psi1_plot_array, Psi2_plot_array = fem_dvr.Plot_Psi_Two_States(wfcnPlot1, wfcnPlot2,plot_title_string="psi_1 and psi_2  with n ="+str(n_Plot),N_plot_points=750,make_plot=True)
    #
    #
    # =============Build Hamiltonian at t=0==============================================
    #     Pass name of potential function explicitly here
    time = 0.0
    H_mat = fem_dvr.Hamiltonian_Two_States(
        pertubation.V_morse_1, pertubation.V_morse_2, Vcoupling, time
    )
    print("\n Completed construction of Hamiltonian at t = 0")
    #
    # ================# Initialize wave function at t = 0 ================================
    #
    #  Bound state of one of the two uncoupled Hamiltonians chosen as initial state
    #
    Cinitial = np.zeros(2 * fem_dvr.nbas, dtype=np.complex)
    for i in range(0, 2 * fem_dvr.nbas):
        Cinitial[i] = (
            wfcnPlot[i] + 0.0 * 1j
        )  # this loop was necessary to make Cinitial an array of complex numbers
    wfcnInitialPlot1 = np.zeros(fem_dvr.nbas, dtype=np.complex)
    wfcnInitialPlot2 = np.zeros(fem_dvr.nbas, dtype=np.complex)
    norm_1 = 0.0
    norm_2 = 0.0
    for j in range(0, fem_dvr.nbas):
        wfcnInitialPlot1[j] = Cinitial[j]
        wfcnInitialPlot2[j] = Cinitial[j + fem_dvr.nbas]
        norm_1 = norm_1 + np.abs(wfcnInitialPlot1[j]) ** 2
        norm_2 = norm_2 + np.abs(wfcnInitialPlot2[j]) ** 2
    #
    # plot initial wave packet
    #
    tinitial = 0.0
    x_Plot_time_array = []
    Psi1_plot_time_array = []
    Psi2_plot_time_array = []
    if want_to_plot_animate is True:
        print("\n Plot initial wave function at t = ", tinitial)
        number_string = str(0.0)
        title = "Wavefunction at t = " + number_string
        (
            x_Plot_array,
            Psi1_plot_array,
            Psi2_plot_array,
        ) = fem_dvr.Plot_Psi_Two_States(
            wfcnInitialPlot1,
            wfcnInitialPlot2,
            plot_title_string=title,
            N_plot_points=750,
            make_plot=True,
        )
        x_Plot_time_array.append(x_Plot_array)
        Psi1_plot_time_array.append(Psi1_plot_array)
        Psi2_plot_time_array.append(Psi2_plot_array)
    #
    # initialize arrays for storing the information for all the frames in the animation of
    # the propagation of the wave function
    times_array = []
    norm_1_array = []
    norm_2_array = []
    times_array.append(tinitial)
    norm_1_array.append(norm_1)
    norm_2_array.append(norm_2)
    #
    # =====================Propagation from t = tinitial to tfinal ============================
    # specify final time and specify number of intervals for plotting
    #
    t_final_propagation = 15.0  # femtoseconds
    tfinal = t_final_propagation / atu_to_fs
    print(
        "\nOverall propagation will be from ",
        tinitial,
        " to ",
        tfinal,
        " = ",
        tfinal * 24.189 / 1000.0,
        "fs in ",
        number_of_time_intervals,
        " intervals",
    )
    #
    #  Loop over n_intervals intervals that make up t = 0 to tfinal
    #
    t_interval = (tfinal - tinitial) / number_of_time_intervals
    time_step = (
        0.1  # atomic time units  -- denser grids require smaller time steps
    )
    for i_time_interval in range(0, number_of_time_intervals):
        t_start = tinitial + i_time_interval * t_interval
        t_finish = tinitial + (i_time_interval + 1) * t_interval
        N_time_steps = np.int(t_interval / time_step)
        Delta_t = t_interval / np.float(N_time_steps)
        #
        # Check norm of initial wavefunction
        #
        norm_initial = 0.0
        norm_1 = 0.0
        for j in range(0, 2 * fem_dvr.nbas):
            norm_initial = norm_initial + np.abs(Cinitial[j]) ** 2
            if j < fem_dvr.nbas:
                norm_1 = norm_1 + np.abs(Cinitial[j]) ** 2
        norm_2 = norm_initial - norm_1
        print(
            "\nNorm of starting  wave function",
            norm_initial,
            " norm on state 1 ",
            norm_1,
            " norm on state 2 ",
            norm_2,
        )
        #
        #  use the propagator from DVR()
        #
        clock_start = timeclock.time()
        if t_start * atu_to_fs < T_pulse:
            Ctfinal = fem_dvr.Crank_Nicolson_Two_States(
                t_start,
                t_finish,
                N_time_steps,
                Cinitial,
                pertubation.V_morse_1,
                pertubation.V_morse_2,
                Vcoupling,
                Is_H_time_dependent=True,
            )
        else:
            Ctfinal = fem_dvr.Crank_Nicolson_Two_States(
                t_start,
                t_finish,
                N_time_steps,
                Cinitial,
                pertubation.V_morse_1,
                pertubation.V_morse_2,
                Vzero,
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
        norm_1 = 0.0
        for j in range(0, 2 * fem_dvr.nbas):
            norm_final = norm_final + np.abs(Ctfinal[j]) ** 2
            if j < fem_dvr.nbas:
                norm_1 = norm_1 + np.abs(Ctfinal[j]) ** 2
        norm_2 = norm_final - norm_1
        print(
            "Norm of final  wave function",
            norm_final,
            " norm on state 1 ",
            norm_1,
            " norm on state 2 ",
            norm_2,
        )
        #
        # Plot of packet at end of each interval using Plot_Psi from DVR()
        #  "False" returns the values but doesn't make a plot, so the values
        #  can be used in the animation.
        #
        print(
            "Plot function propagated to t = ",
            t_finish,
            " atomic time units = ",
            t_finish * 24.189 / 1000.0,
            " fs",
        )
        number_string = str(t_finish)
        title = "Wavefunction at t = " + number_string

        wfcn_at_t_Plot1 = np.zeros(fem_dvr.nbas, dtype=np.complex)
        wfcn_at_t_Plot2 = np.zeros(fem_dvr.nbas, dtype=np.complex)
        for j in range(0, fem_dvr.nbas):
            wfcn_at_t_Plot1[j] = Ctfinal[j]
            wfcn_at_t_Plot2[j] = Ctfinal[j + fem_dvr.nbas]
        # plot using the Plot_Psi_Two_States function in DVRHelper, returns two
        # components of wave function separately
        if want_to_plot_animate is True:
            (
                x_Plot_array,
                Psi1_plot_array,
                Psi2_plot_array,
            ) = fem_dvr.Plot_Psi_Two_States(
                wfcn_at_t_Plot1,
                wfcn_at_t_Plot2,
                plot_title_string=title,
                N_plot_points=750,
                make_plot=False,
            )
            times_array.append(t_finish)
            norm_1_array.append(norm_1)
            norm_2_array.append(norm_2)
            x_Plot_time_array.append(x_Plot_array)
            Psi1_plot_time_array.append(Psi1_plot_array)
            Psi2_plot_time_array.append(Psi2_plot_array)
        #
        #   reset initial packet as final packet for this interval
        for i in range(0, 2 * fem_dvr.nbas):
            Cinitial[i] = Ctfinal[i]
        #
        if want_to_plot_animate is True:
            (
                x_Plot_array,
                Psi1_plot_array,
                Psi2_plot_array,
            ) = fem_dvr.Plot_Psi_Two_States(
                wfcn_at_t_Plot1,
                wfcn_at_t_Plot2,
                plot_title_string=title,
                N_plot_points=750,
                make_plot=True,
            )
            print(
                "Final frame at tfinal = ",
                tfinal,
                "atomic time units  is showing -- ready to prepare animation ",
            )
            print("\n **** close the plot window to proceed ****")
            #
            # print a file with the norms of the components of the wave function as a function of time
            #
            norms_file = open("./Plot_Output/Wavepacket_norms.dat", "w")
            for k in range(0, number_of_time_intervals):
                norms_file.write(
                    "%12.8f %12.8f %12.8f\n"
                    % (times_array[k], norm_1_array[k], norm_2_array[k])
                )
            norms_file.close()

    # ==============================================================================
    # # initialization function: plot the background of each frame
    # ==============================================================================
    def init():
        line.set_data([], [])
        time_text.set_text("")
        ax.set_xlabel(" x (bohr) ", fontsize=16, fontweight="bold")
        ax.set_ylabel(
            " |Psi_1(t)| and |Psi_2(t)| ", fontsize=16, fontweight="bold"
        )
        fig.suptitle(
            "Both Components of Wave Packet", fontsize=16, fontweight="bold"
        )
        ax.plot(
            [x_Plot[0], x_Plot[len(x_Plot) - 1]], [0, 0], "k"
        )  # put in a line at the value Phi = 0
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
        #    re_array = np.real(Psi1_plot_time_array[i])
        #    im_array = np.imag(Psi1_plot_time_array[i])
        abs_array1 = np.abs(Psi1_plot_time_array[i])
        abs_array2 = np.abs(Psi2_plot_time_array[i])
        line1.set_data(x_Plot_array, abs_array1)
        line2.set_data(x_Plot_array, abs_array2)
        #    line3.set_data(x_Plot_array,abs_array)
        time_text.set_text("time = " + time_string + " fs")
        #    return (line1,line2,line3,time_text)
        return (line1, line2, time_text)

    if want_to_plot_animate is True:
        # ==============================================================================
        #  Now plot the animation
        # ==============================================================================
        # reinitialize the figure for the next plot which is the animation
        fig = plt.figure()
        ymax = 2.0
        xmin = x_Plot[0]
        xmax = x_Plot[len(x_Plot) - 1]
        ax = fig.add_subplot(
            111, autoscale_on=False, xlim=(xmin, xmax), ylim=(0, +ymax)
        )
        (line,) = ax.plot([], [], "-r", lw=2)
        (line1,) = ax.plot([], [], "-r", lw=2)
        (line2,) = ax.plot([], [], "-b", lw=2)
        # line3, = ax.plot([], [], 'k',lw=2)
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
        anim.save("Non_adiabatic_coupled_Morse_pots.mp4")
        plt.show()
        print("done")

    if want_to_plot_animate is False:
        print(
            "\n\nSet the command line option want_to_plot_animate=True to see figures, "
            "animation, and create plotting directory.\n\n"
        )


if __name__ == "__main__":
    main()
