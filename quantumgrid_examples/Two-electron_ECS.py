#%%
"""
    C.W. McCurdy 2/2/2021 :bold:`Example Script`

    Time-independent Exterior Complex Scaling (ECS) FEM-DVR example

    Temkin-Poet (s-wave limit) or colinear model of a two-electron atom:
    H- anion or He  bound and autoionizing states

    Uses femdvr.py class library, updated Jan 2021 with 2D routines

    Finite Element Method - Discrete Variable Representation (FEM-DVR)
    for 2D Schrödinger equation using Gauss-Lobatto quadrature in
    each finite element Uses class DVRHelper() to construct FEM-DVR
    points, weights and Kinetic Energy

    The FEM_DVR class  implements Exterior Complex Scaling on
    the FEM-DVR contour.  The value of R0 and the complex scale factor
    :math:`e^{I*theta}` are specified here.  The representation of the potential
    must be able to be evaluated on the complex part of the contour
    or be zero there.

    Example:
        Finds all eigenvalues of complex scaled 2D Hamiltonian
        and plots any one of them, specified by  n_state_plot.
        *NB* Diagonalization does not take advantage of sparse
        banded nature of the 2D Hamiltonian matrix in the FEM-DVR.
        Larger scale practical calculations must do so.

        As an example the :math:`2s^2` autoionizing state of He is chosen
        and the plot of the wave function shows the localized
        resonant state in the middle and autoionization decay
        down the sides parallel to the r_1 and r_2 axes.
        Numerical check: :math:`E_gnd = -2.8790288` for He in s-wave limit
        singlet S autoionizing state :math:`E_res = -0.7228 - 0.001199 i`

        Note that the basis in this example is a simple
        product basis of DVR functions in :math:`r_1` and :math:`r_2`, so
        both singlet (symmetric) and triplet (antisymmetric)
        spatial wave functions appear as eigenfunctions of
        the 2D Hamiltonian.

"""
#%%
# preliminaries to invoke SciPy linear algebra functions
from scipy import linalg as LA

# and NumPy which is used to define pi, sqrt, array, .transpose etc. as np
import numpy as np
import matplotlib.pyplot as plt  # import matplotlib pyplot functions
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
    "--want_to_plot",
    type=click.BOOL,
    default="False",
    help="Set to True if you want to turn on plotting",
)
@click.command()
def main(want_to_plot):

    # Making directories for plots and data
    if want_to_plot is True:
        path = os.getcwd()
        Plot_Output = path + "/Plot_Output_2D"

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
        Data_Output = path + "/Data_Output_2D"
        if os.path.exists(Data_Output):
            print(
                "Directory for output .dat files  already exists", Data_Output
            )
        else:
            print(
                "Attempting to create directory for .dat output files  ",
                Data_Output,
            )
            try:
                os.mkdir(Data_Output)
            except OSError:
                print("Creation of the directory %s failed" % Data_Output)
            else:
                print("Successfully created the directory %s " % Data_Output)
    #%%
    # ===============================FEM_DVR=========================================
    #  Set up the FEM DVR grid given the Finite Element boundaries and order
    #  of Gauss Lobatto quadrature, and compute the Kinetic Energy matrix for
    #  this FEM DVR for Mass = mass set in call (atomic units).
    #
    # Here is where the reduced mass in kinetic energy is set
    Elec_Mass = 1.0
    mu = Elec_Mass
    print("reduced mass mu = ", mu)
    bohr_to_Angstrom = 0.529177
    Hartree_to_eV = 27.211386245988  # NIST ref
    eV_to_wavennumber = 8065.54393734921  # NIST ref on constants + conversions
    Hartree_to_wavenumber = (
        2.1947463136320e5  # value from NIST ref on constants + conversions
    )
    HartreeToKelvin = 315773
    atu_to_fs = 24.18884326509 / 1000
    #
    #  Specify the FEM-DVR grid
    #
    n_order = (
        15  # grid setup here is for calculating bound states of He or H- target
    )
    FEM_boundaries = [0.0, 1.0, 5.0, 10.0, 20.0, 30.0]
    theta_scale = 20
    scale_factor = np.exp(
        1j * theta_scale * np.pi / 180.0
    )  # scale_factor = 1 (Integer) for no scaling
    R0 = 19.0  # Complex scaling begins at greatest FEM boundary < R0 specified here
    print(
        "ECS scaling by ", theta_scale, " degrees at least FEM boundary > ", R0
    )
    fem_dvr = FEM_DVR(
        n_order,
        FEM_boundaries,
        Mass=mu,
        Complex_scale=scale_factor,
        R0_scale=R0,
    )
    print("\nFEM-DVR basis of ", fem_dvr.nbas, " functions")

    # Initialize perturbations
    pertubation = Potential()

    # =============Build one-electron Hamiltonian ===============================
    #     Pass name of one-body potential function explicitly here
    time = 0.0
    H_mat = fem_dvr.Hamiltonian(pertubation.vectorized_V_Coulomb, time)
    print("\n Completed construction of 1D Hamiltonian at t = 0")
    #
    # ========= Calculate two-electron integrals in FEM-DVR basis ===============
    #  These are diagonal in the FEM-DVR basis satisfying
    #     <ij|kl> = delta_ik delta_jl v_ij  See
    #     C W McCurdy et al 2004 J. Phys. B: At. Mol. Opt. Phys. 37 R137
    #
    #  They therefore have the same form as a local potential V(r1,r2)
    #  that would be represented by V(r1_i,r2_j) = v_ij at points on the grid
    #  We pass the matrix v_ij to fem_dvr.Hamiltonian_2D() to build the Hamiltonian
    #  in both cases
    #
    Temkin_Poet = (
        True  # Temkin-Poet model is the radial limit model 1/r12 = 1/r>
    )
    T_mat = 2.0 * fem_dvr.KE_mat
    Two_electron_ints = LA.inv(T_mat)
    R0_ECS = fem_dvr.FEM_boundaries[fem_dvr.i_elem_scale]
    R_max = fem_dvr.FEM_boundaries[len(fem_dvr.FEM_boundaries) - 1]
    if Temkin_Poet:
        print("\n Evaluating two-electron integrals  for Temkin-Poet model")
        print(
            " Grid has R0_ECS = ",
            R0_ECS,
            " R_max ",
            R_max,
            " shape of inverse of KE ",
            Two_electron_ints.shape,
        )
    for i in range(fem_dvr.nbas):
        for j in range(fem_dvr.nbas):
            if Temkin_Poet:
                Two_electron_ints[i, j] = (
                    Two_electron_ints[i, j]
                    / (
                        fem_dvr.x_pts[i + 1]
                        * fem_dvr.x_pts[j + 1]
                        * np.sqrt(fem_dvr.w_pts[i + 1] * fem_dvr.w_pts[j + 1])
                    )
                    + 1.0 / R_max
                )
            else:
                Two_electron_ints[i, j] = pertubation.V_colinear_model(
                    fem_dvr.x_pts[i + 1], fem_dvr.x_pts[j + 1]
                )  # colinear model
    #
    # ============ Build two-electron Hamiltonian ================================
    # pass matrix of two-electron integrals in FEM-DVR or matrix of FEM-DVR
    # representation of local two-electron interaction
    #
    H_mat_2D = fem_dvr.Hamiltonian_2D(H_mat, Two_electron_ints, time)
    print("\n Completed construction of 2D Hamiltonian")
    # =============================================================================
    #
    # ====== Find all the eigenvalues and eigenvectorsof the 2D Hamiltonian =====
    #
    # Here we use LA.eig, which solves the complex generalized eigenvalue & eigenvector problem
    #     (use LA.eigvalsh() for Hermitian matrices)
    #  Called here with arguments right=True for right hand eigenvectors
    #
    print(
        "Calculating ",
        fem_dvr.nbas ** 2,
        " eigenvalues and eigenvectors of 2D Hamiltonian ",
    )
    EigenVals, EigenVecs = LA.eig(
        H_mat_2D, right=True, homogeneous_eigvals=False
    )
    print("after LA.eig()")
    #
    # Sort eigenvalues and eigenvectors using np.argsort() function
    #
    print(" shape of EigenVals on return from LA.eig ", EigenVals.shape)
    print(" shape of EigenVecs on return from LA.eig ", EigenVecs.shape)
    size_eigen_system = EigenVals.shape[
        0
    ]  # note shape is an array two long = (2,dimension)
    print("size_eigen_system = ", size_eigen_system)
    real_parts_eigenvals = np.empty(size_eigen_system)
    for i in range(0, size_eigen_system):
        real_parts_eigenvals[i] = np.real(EigenVals[i])
    idx = np.argsort(
        real_parts_eigenvals
    )  # argsort() returns indices to sort array
    EigenVals = EigenVals[
        idx
    ]  # mysterious python syntax to indirectly address the array
    EigenVecs = EigenVecs[
        :, idx
    ]  # mysterious python syntax to indirectly address 2D array
    #
    # Make a file with sorted complex eigenvalues
    #
    n_energy = fem_dvr.nbas * fem_dvr.nbas
    n_lowest = 10
    print("\n  Lowest ", n_lowest, " eigenvalues of 2D Hamiltonian")
    for i in range(0, min(n_lowest, n_energy)):
        print(
            "E( ",
            i,
            ") =   ",
            np.real(EigenVals[i]),
            "  ",
            np.imag(EigenVals[i]),
            " i",
            "    hartrees",
        )
    #
    if want_to_plot is True:
        file_opened = open(Data_Output + "/Spectrum_ECS.dat", "w")
        for i in range(0, n_energy):
            print(
                np.real(EigenVals[i]),
                "  ",
                np.imag(EigenVals[i]),
                file=file_opened,
            )
        print(" Complete spectrum, sorted on Re(energy) written to file ")
    #
    # =============  Plot of eigenfunctions of 2D Hamiltonian ==================
    #
    n_state_plot = 31  # n_state_plot = 0 is lowest energy (by real part) state
    print(
        "Plotting eigenstate ",
        n_state_plot,
        " with energy ",
        EigenVals[n_state_plot],
    )
    Coef_vector = []
    for i in range(0, n_energy):
        Coef_vector.append(EigenVecs[i, n_state_plot])
    # normalize vector
    sum = 0.0
    for i in range(0, n_energy):
        sum = sum + np.abs(Coef_vector[i]) ** 2
    Norm_const = 1.0 / np.sqrt(sum)
    for i in range(0, n_energy):
        Coef_vector[i] = Norm_const * Coef_vector[i]
    #
    # ========  Optional plot at FEM-DVR points only without interpolation ============
    #
    Plot_at_DVR_only = False
    if want_to_plot is True:
        if Plot_at_DVR_only:
            print(" making plot data file at DVR grid points")
            plot_file_opened = open(Data_Output + "/wave_function_2D.dat", "w")
            for i in range(1, fem_dvr.nbas + 1):
                for j in range(1, fem_dvr.nbas + 1):
                    print(
                        np.real(fem_dvr.x_pts[i]),
                        "   ",
                        np.real(fem_dvr.x_pts[j]),
                        "   ",
                        np.real(
                            Coef_vector[j + (i - 1) * fem_dvr.nbas - 1]
                            / np.sqrt(fem_dvr.w_pts[i] * fem_dvr.w_pts[j])
                        ),
                        file=plot_file_opened,
                    )
                print("    ", file=plot_file_opened)
            exit()
    #
    # =========== Plot of wave function on a grid evenly spaced within FEMs ==============
    #
    if want_to_plot is True:
        plot_pts = 50
        x_Plot_array, y_Plot_array, Psi_plot_array = fem_dvr.Plot_Psi_2D(
            Coef_vector,
            plot_title_string="Real_Psi_on_2D_ECS_grid",
            N_plot_points=plot_pts,
            make_plot=True,
        )

        plot_file_opened = open(Data_Output + "/wave_function_2D.dat", "w")
        print(" Input number of plotting points in each dimension ", plot_pts)
        print(
            "  Actual number of plotting ponts  ",
            len(x_Plot_array),
            " evenly spaced in elements on ECS contour",
        )

        #  write data file for external plotting
        ipt = 0
        plot_pts_on_contour = len(
            x_Plot_array
        )  #  can differ from plot_pts, PLot_Psi_2D puts an integer number in each element
        for i in range(plot_pts_on_contour):
            for j in range(plot_pts_on_contour):
                print(
                    np.real(x_Plot_array[i]),
                    "   ",
                    np.real(y_Plot_array[j]),
                    "   ",
                    np.real(Psi_plot_array[ipt]),
                    file=plot_file_opened,
                )
                ipt = ipt + 1
            print("   ", file=plot_file_opened)
        print(
            " File of points for plot written: ",
            Data_Output + "/wave_function_2D.dat",
        )


if __name__ == "__main__":
    main()
