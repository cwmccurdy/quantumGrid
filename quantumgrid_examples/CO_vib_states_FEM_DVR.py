"""         Chem 210A/B --  C.W. McCurdy -- January 2021

 Read a file with potential curve values and find bound states
              for particle with specified mass

 Finite Element Method - Discrete Variable Representation (FEM-DVR)
 for 1D Schrödinger equation using Gauss-Lobatto quadrature in each finite element
 Uses class FEM_DVR(args) to construct FEM-DVR points, weights and Kinetic Energy

 Shows how to
   (1) Read in and interpolate a potential function known at discrete points
   (2) Use DVRHelper class to build FEM-DVR grid
   (3) Use DVRHelper class to build Hamiltonian in DVR basis
   (4) Find eigenvalues and eigenvectors of Hamiltonian
   (5) Plot eigenfunctions of the Hamiltonian

 Example: CO vibrational states using CI singles, doubles and triples potential curve
          from Psi4
          This potential gives a dissociation energy of :math:`~12.2` eV, not very good
          by comparison to the ~11.1 eV experimental value.
          It yields a :math:`n = 0 -> 1` excitation energy of :math:`2207
          cm^{-1}`
          compared with the value using the NIST values for
          constants of diatomic molecules for :math:`H_2` in the formula
           :math:`E_n = (n+1/2)we - (n+1/2)^2 wexe`, which is :math:`2143
           cm^{-1}`
          So not quite spectroscopic accuracy.

    References
    ----------

    .. bibliography:: _static/refs_examples.bib
      :style: unsrt
"""
# preliminaries to invoke SciPy linear algebra functions
from scipy import linalg as LA

# and NumPy which is used to define pi, sqrt, array, .transpose etc. as np
import numpy as np
import matplotlib.pyplot as plt  # import matplotlib pyplot functions
import os  # functions to manipulate files and directories

# Needed to read in data files distributed in quantumgrid examples
from pathlib import Path

# for debugging
# import sys
# sys.path.append("../")

# Import our module!
from quantumgrid.femdvr import FEM_DVR
from quantumgrid.potential import Potential

import click


@click.option(
    "--want_to_plot",
    type=click.BOOL,
    default="False",
    help="Set to True if you want to turn on plotting",
)
@click.command()
def main(want_to_plot):
    #
    # ============== Make Directory for Plots if it's not there already ==========
    #
    # detect the current working directory and print it
    path = os.getcwd()
    print("The current working directory is %s" % path)
    # define the name of the directory to be created
    Plot_Output = path + "/Plot_Output"

    if want_to_plot is True:
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
    #
    # === Make Directory for for the .dat output  if it's not there already ======
    #
    # detect the current working directory and print it
    path = os.getcwd()
    print("The current working directory is %s" % path)
    # define the name of the directory to be created
    Data_Output = path + "/Data_Output"
    if os.path.exists(Data_Output):
        print("Directory for output .dat files  already exists", Plot_Output)
    else:
        print(
            "Attempting to create directory for .dat output files  ",
            Plot_Output,
        )
        try:
            os.mkdir(Data_Output)
        except OSError:
            print("Creation of the directory %s failed" % Data_Output)
        else:
            print("Successfully created the directory %s " % Data_Output)
    # ===================================================================
    if want_to_plot is True:
        wfcn_plotfile = open(
            "Data_Output/wavefunctions.dat", "w"
        )  # data file for saving wavefunctions
    #
    # =============Constants and conversion factors ==============
    Daltons_to_eMass = 1822.89
    bohr_to_Angstrom = 0.529177
    Hartree_to_eV = 27.211386245988  # NIST ref
    eV_to_wavennumber = 8065.54393734921  # NIST ref on constants + conversions
    Hartree_to_wavenumber = (
        2.1947463136320e5  # value from NIST ref on constants + conversions
    )
    atu_to_fs = 24.18884326509 / 1000
    HartreeToKelvin = 315773
    #
    # =====================================FEM_DVR===================================
    #  Set up the FEM DVR grid given only the Finite Element boundaries and order
    #  of Gauss Lobatto quadrature,  and compute the Kinetic Energy matrix for
    #  this FEM DVR for Mass = mass set in call (atomic units).
    #
    # Here is where the reduced mass (or mass for any 1D problem) is set
    O_Mass = 15.99491461957  # O 16 mass from NIST tables
    C_Mass = 12.0  # C 12 mass
    mu = (O_Mass * C_Mass / (O_Mass + C_Mass)) * Daltons_to_eMass
    n_order = 30
    FEM_boundaries = [
        0.801,
        1.25,
        1.5,
        2.0,
        2.5,
        3.0,
        3.5,
        4.0,
        4.5,
        5.0,
        5.5,
        6.0,
        6.5,
        7.0,
        7.5,
        8.0,
        8.5,
        8.99,
    ]
    fem_dvr = FEM_DVR(n_order, FEM_boundaries, Mass=mu)
    print("\nFEM-DVR basis of ", fem_dvr.nbas, " functions")
    #
    #   Function to define potential at x and t (if potential is time-dependent)
    #   DVRHelper class library expects function for V(x,t) in general
    #
    # =================================Potential=====================================
    #  Read in files with points for potential curve in hartrees
    #  and load in arrays for interpolation
    path = Path(__file__).parent.absolute()
    perturbation = Potential(path / "potcurve_CISDT_CO_ccpvDZ.dat")

    n_vals_pot = perturbation.r_data.shape[0]
    #

    # ===========================================================================================
    #  Plot potential on the DVR grid points on which the wavefunction is defined
    #  and ALSO the interpolation to check we are using the potential that we mean to.
    #
    if want_to_plot is True:
        print("\n Plot potential ")
        x_Plot = []
        pot_Plot = []
        for j in range(0, n_vals_pot):
            x_Plot.append(np.real(perturbation.r_data[j]))
            pot_Plot.append(np.real(perturbation.V_data[j]))
        plt.suptitle("V(r) interpolation", fontsize=14, fontweight="bold")
        string = "V input points"
        plt.plot(x_Plot, pot_Plot, "ro", label=string)
        #
        x_Plot = []
        pot_Plot = []
        Number_plot_points = 731
        dx = (
            np.real(fem_dvr.x_pts[fem_dvr.nbas - 1]) - np.real(fem_dvr.x_pts[0])
        ) / float(Number_plot_points - 1)
        time = 0.0  # dummy time in general call to potential function
        for j in range(0, Number_plot_points):
            x = np.real(fem_dvr.x_pts[0]) + j * dx
            try:
                x >= perturbation.r_data[0] and x <= perturbation.r_data[
                    n_vals_pot - 1
                ]
            except IndexError:
                print(
                    "Number of plot points is out of range of pertubation data"
                )
            x_Plot.append(x)
            pot_Plot.append(perturbation.V_Interpolated(x, time))

        plt.plot(
            x_Plot, pot_Plot, "-b", label="Interpolation on DVR grid range"
        )
        plt.legend(loc="best")
        plt.xlabel(" x ", fontsize=14)
        plt.ylabel("V", fontsize=14)
        print(
            "\n Running from terminal, close figure window to proceed and make .pdf file of figure"
        )
        #   Insert limits if necessary
        # xmax = float(rmax)  # CWM: need to use float() to get plt.xlim to work to set x limits
        # plt.xlim([0,xmax])
        # number_string = str(a)
        plt.savefig(
            "Plot_Output/" + "Plot_potential" + ".pdf", transparent=False
        )
        plt.show()
    #
    #
    # =============Build Hamiltonian (using general routine with dummy time t=0)=========
    #     Pass name of potential function explicitly here
    time = 0.0
    H_mat = fem_dvr.Hamiltonian(perturbation.vectorized_V_Interpolated, time)
    print("\n Completed construction of Hamiltonian ")
    # ====================================================================================
    #
    # Find all the eigenvalues of the Hamiltonian so we can compare with known bound state energies
    # or make a plot of the spectrum -- For a time-independent Hamiltonian example here
    #
    EigenVals = LA.eigvalsh(H_mat)
    #
    n_energy = 20
    for i in range(0, n_energy):
        print(
            "E( ",
            i,
            ") =   ",
            EigenVals[i],
            " hartrees, excitation energy = ",
            (EigenVals[i] - EigenVals[0]) * Hartree_to_wavenumber,
            " cm^-1",
        )

    if want_to_plot is True:
        # ====================================================================================
        #
        # Extract the n_Plot'th eigenfunction for plotting and use as initial wave function
        # to test propagation
        #
        number_of_eigenvectors = 20
        #
        #  Here n_Plot picks which eigenfunction to plot
        n_Plot = 10  # pick a state of this potential to plot < number_of_eigenvectors -1
        #
        print(
            "Calculating ",
            number_of_eigenvectors,
            " eigenvectors for plotting eigenfunctions",
        )
        EigenVals, EigenVecs = LA.eigh(
            H_mat, eigvals=(0, number_of_eigenvectors)
        )
        wfcnPlot = []
        for j in range(0, fem_dvr.nbas):
            wfcnPlot.append(EigenVecs[j, n_Plot])
        #
        # normalize  wave function from diagonalization
        #
        norm_squared = 0.0
        for j in range(0, fem_dvr.nbas):
            norm_squared = norm_squared + np.abs(wfcnPlot[j]) ** 2
        wfcnPlot = wfcnPlot / np.sqrt(norm_squared)
        norm_squared = 0.0
        for j in range(0, fem_dvr.nbas):
            norm_squared = norm_squared + np.abs(wfcnPlot[j]) ** 2
        print("Norm of wave function being plotted is ", np.sqrt(norm_squared))
        #
        # ================# Plot the  wave function specified by n_Plot above======================
        #   It must be type np.complex for this general wave function plotting logic
        #
        Cinitial = np.zeros((fem_dvr.nbas), dtype=np.complex)
        wfcnInitialPlot = np.zeros((fem_dvr.nbas), dtype=np.complex)
        for j in range(0, fem_dvr.nbas):
            Cinitial[j] = wfcnPlot[j]
        #
        print(
            "\n Plot wave function ",
            n_Plot,
            " (numbered in order of increasing energy)",
        )
        title = "Wavefunction" + str(n_Plot)
        #  note that the dvr.Plot_Psi function makes a .pdf file in the Plot_Output directory
        #  That's what make_plot=True controls.
        x_Plot_array, Psi_plot_array = fem_dvr.Plot_Psi(
            Cinitial,
            plot_title_string=title,
            N_plot_points=Number_plot_points,
            make_plot=True,
        )
        # write the data in file also
        for j in range(0, Number_plot_points):
            print(
                x_Plot_array[j],
                "  ",
                np.real(Psi_plot_array[j]),
                "  ",
                np.imag(Psi_plot_array[j]),
                file=wfcn_plotfile,
            )
        print("&  \n ", file=wfcn_plotfile)
        #
        # plot square of wave function (radial probability distribution)
        #
        Csquared = np.zeros(fem_dvr.nbas, dtype=np.complex)
        raverage = 0.0
        for i in range(fem_dvr.nbas):
            Csquared[i] = (Cinitial[i] ** 2) / np.sqrt(fem_dvr.w_pts[i + 1])
            #   compute <r> for this wave function
            #   note that Cinitial[i] contains the necessary weight, sqrt(dvr.w_pts[i+1])
            raverage = raverage + Cinitial[i] ** 2 * fem_dvr.x_pts[i + 1]
        print(
            "\n Average value of r using DVR for the integral, <r> = ", raverage
        )
        title = "Wavefunction" + str(n_Plot) + "^2"
        #  note that the dvr.Plot_Psi function makes a .pdf file in the Plot_Output directory
        #  That's what make_plot=True controls.
        x_Plot_array, Psi_plot_array = fem_dvr.Plot_Psi(
            Csquared,
            plot_title_string=title,
            N_plot_points=Number_plot_points,
            make_plot=True,
        )
        # write the data in file also
        for j in range(0, Number_plot_points):
            print(
                x_Plot_array[j],
                "  ",
                np.real(Psi_plot_array[j]),
                "  ",
                np.imag(Psi_plot_array[j]),
                file=wfcn_plotfile,
            )
        print("&  \n ", file=wfcn_plotfile)
        #


if __name__ == "__main__":
    main()
