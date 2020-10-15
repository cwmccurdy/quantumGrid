"""
    C.W. McCurdy 03/13/2020
    Time-independent Exterior Complex Scaling (ECS) FEM-DVR example
    Uses femdvr.py class and potential class.

    Finite Element Method - Discrete Variable Representation (FEM-DVR)
    for 1D Schroedinger equation using Gauss-Lobatto quadrature in
    each finite element Uses class DVRHelper() to construct FEM-DVR
    points, weights and Kinetic Energy

    For time-independent potential, this example implements Exterior
    Complex Scaling on the FEM-DVR contour.  The value of R0 and the
    complex scale factor :math:`e^{I*theta}` are specified.  The representation
    of the pootential must be able to be evaluated on the complex part
    of the contour.

    Example:
       Finds all eigenvalues of complex scaled Hamiltonian and
       plots any one of them, specified by n_Plot

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

# and NumPy which is used to define pi, sqrt, array, .transpose etc. as
import numpy as np
import matplotlib.pyplot as plt  # import matplotlib pyplot functions
from matplotlib import animation  # for animation from same class library
import os  # functions to manipulate files and directories
import time as timeclock  # for timing parts of the calculation during debugging

# Needed to import our classes
# import sys

# sys.path.append("../")

# Importing our classes
from quantumgrid.femdvr import FEM_DVR
from quantumgrid.potential import Potential

import sys
import click


@click.command()
def main(args=None, function=None):
    #
    # ================ Make Directory for Plots if it's not there already =============
    #
    # detect the current working directory and print it
    path = os.getcwd()
    print("The current working directory is %s" % path)
    # define the name of the directory to be created
    Plot_Output = path + "/Plot_Output"
    if os.path.exists(Plot_Output):
        print("Directory for wave function plots already exists", Plot_Output)
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
    n_order = 25
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
        15.0,
        20.0,
        25.0,
        30.0,
        50,
    ]
    # parameters for Fig 2 in Turner, McCurdy paper
    # Julia Turner and C. William McCurdy, Chemical Physics 71(1982) 127-133
    scale_factor = np.exp(1j * 37.0 * np.pi / 180.0)
    R0 = 22.75
    fem_dvr = FEM_DVR(
        n_order,
        FEM_boundaries,
        Mass=mu,
        Complex_scale=scale_factor,
        R0_scale=R0,
    )
    print("\nFEM-DVR basis of ", fem_dvr.nbas, " functions")
    #

    pertubation = Potential()
    # ==================================================================================
    #  Plot potential on the DVR grid points on which the wavefunction is defined
    print("\n Plot potential ")
    print("Test V", pertubation.V_Bernstein(5.0 + 0.0 * 1j, 0.0))
    print("Test V", pertubation.V_Bernstein(10.0 + 0.5 * 1j, 0.0))
    time = 0.0
    x_Plot = []
    pot_Plot = []
    for j in range(0, fem_dvr.nbas):
        x_Plot.append(np.real(fem_dvr.x_pts[j + 1]))
        pot_Plot.append(
            np.real(pertubation.V_Bernstein(fem_dvr.x_pts[j + 1], time))
        )
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
    plt.savefig("Plot_Output/" + "Plot_potential" + ".pdf", transparent=False)
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
    n_Plot = 292
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
    # exit()

    # ====================================================================================
    #
    # Extract the n_Plot'th eigenfunction for plotting
    #
    # pick one of the eigenstates of Potential to plot
    # numbering can depend on numpy and python installation that determines
    # behavior of the linear algebra routines.
    n_Plot = 292
    print(
        "Calculating ",
        fem_dvr.nbas,
        " eigenvectors for plotting eigenfunctions",
    )
    EigenVals2, EigenVecs = LA.eig(H_mat, right=True, homogeneous_eigvals=True)
    wfcnPlot = []
    for j in range(0, fem_dvr.nbas):
        wfcnPlot.append(EigenVecs[j, n_Plot])
    #
    # normalize  wave function from diagonalization
    # using integration of square on the contour
    # following the old Rescigno McCurdy idea for partial widths
    norm_squared = 0.0
    for j in range(0, fem_dvr.nbas):
        norm_squared = norm_squared + (wfcnPlot[j]) ** 2
    wfcnPlot = wfcnPlot / np.sqrt(norm_squared)
    norm_squared = 0.0
    for j in range(0, fem_dvr.nbas):
        norm_squared = norm_squared + (wfcnPlot[j]) ** 2
    print(
        "Norm of wave function from int psi^2 on contour being plotted is ",
        np.sqrt(norm_squared),
    )
    # Plot wave function -- It must be type np.complex
    Cinitial = np.zeros((fem_dvr.nbas), dtype=np.complex)
    wfcnInitialPlot = np.zeros((fem_dvr.nbas), dtype=np.complex)
    for j in range(0, fem_dvr.nbas):
        Cinitial[j] = wfcnPlot[j]
    #
    # plot n_Plot'th eigenfunction
    #
    print(
        "\n Plot Hamiltonian eigenfunction number ",
        n_Plot,
        " with energy ",
        EigenVals[0, n_Plot],
    )
    number_string = str(n_Plot)
    title = "Wavefunction number = " + number_string
    wfn_plot_points = 1000
    x_Plot_array, Psi_plot_array = fem_dvr.Plot_Psi(
        Cinitial,
        plot_title_string=title,
        N_plot_points=wfn_plot_points,
        make_plot=True,
    )
    #
    filename = "wavefunction" + number_string + ".dat"
    file_opened = open(filename, "w")
    print_points = len(x_Plot_array)
    # print("x_Plot_array shape ",print_points)
    for i in range(print_points):
        print(
            np.real(x_Plot_array[i]),
            "  ",
            np.imag(x_Plot_array[i]),
            "  ",
            np.real(Psi_plot_array[i]),
            "  ",
            np.imag(Psi_plot_array[i]),
            file=file_opened,
        )


if __name__ == "__main__":
    sys.exit(main())  # pragma: no covers
