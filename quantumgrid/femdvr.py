"""Finite Eement Method -- Discrete Variable Representation (DVR) for
   1D Schroedinger equation using Gauss-Lobatto quadrature
   C. William McCurdy, Zachary Streeter, and Giuseppe Barbalinardo -- UC Davis

   Includes time propagation with Crank-Nicolson propagator
   and routines for dynamics on two coupled potential curves.

   November 2019
    time-dependent potentials allowed for two-state
    propagation vectorized logic in Hamiltonian and Potential builds

   March 2020
    Exterior Complex Scaling implemented for single state portion
    both time-independent and time-dependent calculations
    Vectorized logic for diagonals of Hamiltonian can be commented and
    alternate logic uncommented for potential functions that don't
    vectorize correctly.  2-state routines are not yet implemented for ECS

"""
import time as timeclock  # for timing parts of the calculation

import matplotlib.pyplot as plt  # for the wave function plotting function

# Import  NumPy which is used to define pi, sqrt, array, .transpose etc. as
import numpy as np
from numpy import linalg as LA  # for linear eq solve in Crank-Nicolson propagator

#                           Preliminaries
#
# import the gauss_lobatto routine from sympy -- basis of DVR
from sympy.integrals.quadrature import gauss_lobatto


class FEM_DVR(object):
    def __init__(self, n_order, FEM_boundaries, Mass=1, Complex_scale=1, R0_scale=0.0):
        """Constructor method
        Builds 1D  FEM DVR grid and kinetic energy matrix representation in
        normalized FEM DVR basis of shape functions and bridging functions
        Finite elements of any sizes determined by FEM_boundaries array, but
        all with same order DVR

        Args:
            n_order (int): The DVR order. More dense takes longer but improves
                            acuracy
            FEM_boundaries (ndarray): An array of 'double' or 'complex' boundaries for the
                                        Finite-Elements.
            Complex_scale (int,optional): Flag that starts complex scaling at grid
                                    boundary closest to R0 if .ne. 1. Defaults
                                    to 1
            R0_scale (complex,optional): Grid point where the complex tail starts.
                                    Defaults to 0.0

        Attributes:
            n_order (int): Can access the DVR order later
            FEM_boundaries (ndarray): Can access the 'double' or 'complex'
                Finite-Element boundaries later
            nbas (int): Total number of basis functions in the FEM-DVR grid
            x_pts (ndarray): Array of 'double' or 'complex' DVR points
            w_pts (ndarray): Array of 'double' or 'complex' DVR weight
            KE_mat (ndarray): Kintetic Energy matrix, dimension is
                nbas x nbas of 'double' or 'complex' types
            i_elem_scale (int): index of last real DVR point

        """
        self.n_order = n_order
        self.FEM_boundaries = FEM_boundaries
        N_elements = len(FEM_boundaries) - 1
        print(
            "\nFrom FEM_DVR_build: building grid and KE with Gauss Lobatto quadrature of order ",
            n_order,
        )
        #
        # Complex the FEM_boundaries if Complex_scale .ne. 1
        #
        i_elem_scale = 0.0
        if Complex_scale != 1:
            for i_elem in range(0, N_elements):
                if FEM_boundaries[i_elem] >= R0_scale:
                    R0 = FEM_boundaries[i_elem]
                    i_elem_scale = i_elem
                    break
            for i_elem in range(i_elem_scale, N_elements + 1):
                FEM_boundaries[i_elem] = R0 + Complex_scale * (
                    FEM_boundaries[i_elem] - R0
                )
        #
        for i_elem in range(0, N_elements):
            print(
                "element ",
                i_elem + 1,
                " xmin = ",
                FEM_boundaries[i_elem],
                " xmax = ",
                FEM_boundaries[i_elem + 1],
            )
        #
        #  sympy gauss_lobatto(n,p) does extended precision, apparently symbolically, to
        #  produce points and weights on (-1,+1).  It is slow for higher  orders
        n_precision = 18
        xlob, wlob = gauss_lobatto(n_order, n_precision)
        #
        w_lobatto = []
        x_lobatto = []
        # make the elements of the extended precision points and weights from sympy
        # into ordinary floating point numbers for subsequent use
        for i in range(0, n_order):
            w_lobatto.append(float(wlob[i]))
            x_lobatto.append(float(xlob[i]))
        # create temporary array to hold kinetic energy before first and last points of grid are removed
        # to enforce boundary condition psi = 0 at ends of grid
        # all initializations with np.zeros are complex for ECS
        full_grid_size = n_order * N_elements - (N_elements - 1)
        KE_temp = np.zeros((full_grid_size, full_grid_size), dtype=np.complex)
        # Build FEM_DVR and Kinetic Energy matrix in same for loop:
        x_pts = np.zeros(full_grid_size, dtype=np.complex)
        w_pts = np.zeros(full_grid_size, dtype=np.complex)
        for i_elem in range(0, N_elements):
            xmin = FEM_boundaries[i_elem]
            xmax = FEM_boundaries[i_elem + 1]
            x_elem = []
            w_elem = []
            for i in range(0, n_order):
                x_elem.append(((xmax - xmin) * x_lobatto[i] + (xmax + xmin)) / 2.0)
                w_elem.append((xmax - xmin) * w_lobatto[i] / 2.0)
            for i in range(0, n_order):
                shift = i_elem * n_order
                if i_elem > 0:
                    shift = i_elem * n_order - i_elem
                x_pts[i + shift] = x_elem[i]
                w_pts[i + shift] = w_pts[i + shift] + w_elem[i]
            # DVR().Kinetic_energy_FEM_block gives -1/2 d^x/dx^2 in unnormalized DVR basis
            KE_elem = self.Kinetic_Energy_FEM_block(n_order, x_elem, w_elem)
            for i in range(0, n_order):
                for j in range(0, n_order):
                    KE_temp[
                        i + (i_elem * (n_order - 1)), j + (i_elem * (n_order - 1)),
                    ] = (
                        KE_temp[
                            i + (i_elem * (n_order - 1)), j + (i_elem * (n_order - 1)),
                        ]
                        + KE_elem[i, j]
                    )
        nbas = full_grid_size - 2
        KE_mat = np.zeros((nbas, nbas), dtype=np.complex)
        # apply normalizations of DVR basis to KE matrix, bridging fcns have w_left + w_right
        for i in range(0, nbas):
            for j in range(0, nbas):
                KE_mat[i, j] = KE_temp[i + 1, j + 1] / np.sqrt(
                    w_pts[i + 1] * w_pts[j + 1]
                )
        # KE matrix complete, return full grid with endpoints and KE_mat which is
        # nbas x nbas
        self.nbas = nbas
        self.x_pts = x_pts
        self.w_pts = w_pts
        self.KE_mat = KE_mat
        self.rescale_kinetic(Mass)
        self.i_elem_scale = i_elem_scale

    def deriv(self, n, x, w):
        """
        Derivative of Unormalized lobatto shape functions. If FEM-DVR grid has
        complex tail, then derivative array is complex

        Args:
            n (int): the dimension of the derivative matrix
            x (ndarray): array of 'double' or 'complex' DVR points
            w (ndarray): array of 'double' or 'complex' DVR weights

        Returns:
            deriv_matrix (ndarray): square matrix of size n x n
        """
        deriv_matrix = np.identity(n, dtype=np.complex)
        deriv_matrix = 0 * deriv_matrix
        for i in range(0, n):
            der = 0.0
            for j in range(0, n):
                if j != i:
                    der = der + 1.0 / (x[i] - x[j])
            deriv_matrix[i, i] = der
        for i in range(0, n - 1):
            for j in range(i + 1, n):
                de = 1.0 / (x[i] - x[j])
                for k in range(0, n):
                    if k != i and k != j:
                        de = de * (x[j] - x[k]) / (x[i] - x[k])
                deriv_matrix[i, j] = de
                deriv_matrix[j, i] = -deriv_matrix[i, j] * w[j] / w[i]
        #  Uncomment this logic to change function to
        #  scale to give result for normalized DVR basis fcns
        #       for i in range(0,n):
        #          term = 1.0/np.sqrt(w[i])
        #          for j in range(0,n):
        #              deriv_matrix[i,j] = deriv_matrix[i,j]*term
        return deriv_matrix

    def Kinetic_Energy_FEM_block(self, n, x, w):
        """
        Calculates block of Kinetic energy from one Finite-Element. If FEM-DVR grid has
        complex tail, then Kmat is complex

        Args:
            n (int): the dimension of the derivative matrix
            x (ndarray): array of 'double' or 'complex' DVR points
            w (ndarray): array of 'double' or 'complex' DVR weights

        Returns:
            Kmat (ndarray): square matrix of size n x n
        """
        Dmat = self.deriv(n, x, w)
        Kmat = np.identity(n, dtype=np.complex)
        Kmat = 0.0 * Kmat
        for i in range(0, n):
            for j in range(0, n):
                sum = 0.0
                for k in range(0, n):
                    sum = sum + Dmat[i, k] * Dmat[j, k] * w[k]
                Kmat[i, j] = sum / 2.0
        return Kmat

    def Hamiltonian(self, V_potential, time):
        """
        Build Hamiltonian given Kinetic energy matrix and potential function
        Add potential to diagonal -- specified by V_potential(x,t)

        Args:
            V_potential (function): A potential function the caller must provide
            time (int): Time dependence of the Hamiltonian. This argument is
                passed to the potential to simulate turning on a field pertubation,
                for example

        Returns:
            H_mat (ndarray): Hamiltonian for the system defined by the caller provided V_potential, size nbas x nbas
        """
        nbas = self.nbas
        KE_mat = self.KE_mat
        x = self.x_pts
        H_mat = KE_mat.copy()
        # this vectorized logic can be  replaced
        # j = np.arange(nbas)
        # H_mat[j, j] = H_mat[j, j] + V_potential(x[j + 1], time)  #  Potential added on diagonal
        #
        # By this logic for potential functions V_potential() that don't vectorize properly
        #
        for j in range(nbas):
            # Potential added on diagonal
            H_mat[j, j] = H_mat[j, j] + V_potential(np.real(x[j + 1]), time)
        return H_mat

    def Potential_Two_States(self, V_potential_1, V_potential_2, V_coupling, time):
        """
        Build potential function

        Diagonal blocks are FEM-DVR Hamiltonians for V_potential_1 and Vpotential_2
        Off diagonal blocks are diagonal representation of V_coupling(x,t)
        Potential functions are passed in

        Args:
            V_potential_1 (function): A potential function for the first state, the
                caller must provide
            V_potential_2 (function): A potential function for the second state,
                the caller must provide
            V_coupling (function): A potential coupling function providing the
                off-diagnol blocks of the returned potential,  the caller must provide
            time (int): Time dependence of the Hamiltonian. This argument is
                passed to the potential to simulate turning on a field pertubation,
                for example

        Returns:
            potential (ndarray): potential for the two state system defined by the caller provided V_potential_1, V_potential_2, and V_coupling, size nbas x nbas
        """
        nbas = self.nbas
        potential = np.zeros((2 * nbas, 2 * nbas), dtype=np.complex)
        x = self.x_pts
        potential_1 = np.zeros((nbas, nbas))
        potential_2 = np.zeros((nbas, nbas))

        jj = np.arange(nbas)
        potential_1[jj, jj] = V_potential_1(x[jj + 1], time)  # potential for state 1
        potential_2[jj, jj] = V_potential_2(x[jj + 1], time)  # potential for state 2
        #  Load potentials into full potential matrix, vectorized logic
        #  The potentials are placed  on the diagonals and the coupling on the
        #   diagonals of the coupling blocks
        ii = np.arange(nbas)
        potential[ii, ii] = potential_1[ii, ii]
        potential[ii + nbas, ii + nbas] = potential_2[ii, ii]
        kk = np.arange(nbas)
        potential[kk, kk + nbas] = V_coupling(x[kk + 1], time)
        potential[kk + nbas, kk] = np.conj(V_coupling(x[kk + 1], time))
        return potential

    def Hamiltonian_Two_States(self, V_potential_1, V_potential_2, V_coupling, time):
        """
        Build Hamiltonian given Kinetic energy matrix and potential function

        Diagonal blocks are FEM-DVR Hamiltonians for V_potential_1 and Vpotential_2
        Off diagonal blocks are diagonal representation ov V_coupling(x,t)
        Potential functions are passed in

        Args:
            V_potential_1 (function): A potential function for the first state, the
                caller must provide
            V_potential_2 (function): A potential function for the second state,
                the caller must provide
            V_coupling (function): A potential coupling function providing the
                off-diagnol blocks of the returned potential,  the caller must provide
            time (int): Time dependence of the Hamiltonian. This argument is
                passed to the potential to simulate turning on a field pertubation,
                for example

        Returns:
            H_mat (ndarray): Hamiltonian for the two state system defined by the caller provided V_potential_1, V_potential_2, and V_coupling, size nbas x nbas
        """
        #
        #  Build Hamiltonian given Kinetic energy matrix and potential function
        #
        #  Diagonal blocks are FEM-DVR Hamiltonians for V_potential_1 and Vpotential_2
        #  Off diagonal blocks are diagonal representation ov V_coupling(x,t)
        #  Potential functions are passed in
        #
        nbas = self.nbas
        KE_mat = self.KE_mat
        H_mat = self.Potential_Two_States(
            V_potential_1, V_potential_2, V_coupling, time
        )
        # Vectorized logic to add kinetic energy blocks to the diagonals
        j = np.arange(nbas)
        H_mat[0:nbas, j] += KE_mat[0:nbas, j]
        H_mat[nbas : 2 * nbas, j + nbas] += KE_mat[0:nbas, j]
        return H_mat

    def psi_evaluate(self, x, coef_vector):
        """
        Evaluate a function represented by a vector of coefficients of the
        FEM-DVR basis functions at the point x.  Array containing vector of
        coefficients does not contain coefficients for the beginning and end
        of the FEM-DVR grid where boundary condx. enforce wavefunction =  zero.

        Find which finite element x is in (or on the boundary of)

        Args:
            x (complex): FEM-DVR grid point of evalutation of :math:`\Psi`
            coef_vector (ndarray): Function representation in the FEM-DVR basis

        Returns:
            psi_value (complex): The value of :math:`\Psi(x)`
        """
        x_grid = self.x_pts
        w_grid = self.w_pts
        FEM_boundaries = self.FEM_boundaries
        N_order = self.n_order

        N_elements = len(FEM_boundaries) - 1
        i_elem = -1
        for i in range(0, N_elements):
            if (
                np.real(x) >= np.real(FEM_boundaries[i])
                and np.real(x) <= np.real(FEM_boundaries[i + 1]) + 1.0e-9
            ):
                i_elem = i + 1
        if i_elem < 0:
            print(
                "x value ",
                x,
                " is out of range ",
                FEM_boundaries[0],
                ", ",
                FEM_boundaries[N_elements],
                " in psi_evaluate()",
            )
            exit()
        # elements are numbered 1 through N_elements
        # beginning and ending indices in x_grid and w_grid of this element
        index_left = (i_elem - 1) * N_order - (i_elem - 1)
        #  loop on lobatto shape functions in this interval including endpoints
        #  evaluate each shape function at x and multiply by coefficient to accumulate
        #  sum over all contributions.
        sum_val = 0.0
        for j_fcn in range(index_left, index_left + N_order):
            if j_fcn == 0 or j_fcn == len(x_grid) - 1:
                Coef = 0.0  # psi is assumed to be zero on the endpoints of the grid
            else:
                Coef = coef_vector[j_fcn - 1]
            funcval = 1 / np.sqrt(w_grid[j_fcn])
            for i in range(index_left, index_left + N_order):
                if i != j_fcn:
                    funcval = funcval * (x - x_grid[i]) / (x_grid[j_fcn] - x_grid[i])
            sum_val = sum_val + Coef * funcval
        psi_value = sum_val
        return psi_value

    def Plot_Psi(
        self,
        Psi_coefficients,
        plot_title_string="Plot of FEM-DVR representation",
        N_plot_points=500,
        make_plot=True,
    ):
        """
        Quick generic plot of function represented by FEM DVR coefficients

        Args:
            Psi_coefficients (ndarray): Array of type 'double' or 'complex'
                coefficients for the representation of :math:`\Psi` on the FEM-DVR grid
            plot_title_string (string): Title of plot, defaults to "Plot of
                FEM-DVR representation"
            N_plot_points (int): Number of points to plot, default 500
            make_plot (bool): Boolean that turns off/on this plotting feature,
                default true

        Returns:
            x_Plot, Psi_plot (ndarray, ndarrya): x and y coordinates of the graph
        """
        #
        # quick generic plot of function represented by FEM DVR coefficients
        # default number of points to plot specified in keyword argument
        # returns x, y points for a plot without making it if make_plot=False
        x_grid = self.x_pts
        w_grid = self.w_pts
        FEM_boundaries = self.FEM_boundaries
        N_order = self.n_order
        N_elements = len(FEM_boundaries) - 1
        # dx = (FEM_boundaries[N_elements] - FEM_boundaries[0]) / (N_plot_points - 1)
        # print("dx step in Plot_Psi  ",dx)
        Psi_Plot = []
        x_Plot = []
        # Build array of plot points on the contour in x, ECS contour if complex scaling is on
        N_pts_per_elem = np.int(N_plot_points / N_elements)
        for i_elem in range(0, N_elements):
            dx = (FEM_boundaries[i_elem + 1] - FEM_boundaries[i_elem]) / N_pts_per_elem
            if i_elem == N_elements - 1:
                dx = (FEM_boundaries[i_elem + 1] - FEM_boundaries[i_elem]) / (
                    N_pts_per_elem - 1
                )
            for k_pt in range(0, N_pts_per_elem):
                xval = FEM_boundaries[i_elem] + k_pt * dx
                x_Plot.append(xval)
                Psi_Plot.append(self.psi_evaluate(xval, Psi_coefficients))
        #  Points and values of Psi are calculated, now make the plot
        if make_plot:
            figure = plt.figure()
            plt.suptitle(plot_title_string, fontsize=14, fontweight="bold")
            string1 = "Re(psi(t))"
            string2 = "Im(psi(t))"
            string3 = "Abs(psi(t))"
            plt.plot(np.real(x_Plot), np.real(Psi_Plot), "-r", label=string1)
            plt.plot(np.real(x_Plot), np.imag(Psi_Plot), "-g", label=string2)
            plt.plot(np.real(x_Plot), np.abs(Psi_Plot), "-k", label=string3)
            plt.legend(loc="best")
            plt.xlabel(" x ", fontsize=14)
            plt.ylabel("psi", fontsize=14)
            #     limits if necessary
            # xmax = float(rmax)  # CWM: need to use float() to get plt.xlim to work to set x limits
            # plt.xlim([0,xmax])
            # rmax=16.0
            # x_max_plot = float(rmax)  # CWM: need to use float() to get plt.xlim to work to set x limits
            # plt.xlim([0,x_max_plot])
            print(
                "\n Running from terminal, close figure window to proceed and make .pdf file of figure"
            )
            plt.savefig("Plot_Output/" + plot_title_string + ".pdf", transparent=False)
            plt.show()  # note plt.show() evidently clears everything for this plot
        return x_Plot, Psi_Plot  # returns the x, y coordinates for a graph

    def Plot_Psi_Two_States(
        self,
        Psi_coefficients1,
        Psi_coefficients2,
        plot_title_string="Plot of FEM-DVR representation",
        N_plot_points=500,
        make_plot=True,
    ):
        """
        Quick generic plot of wave function represented by FEM DVR coefficients
        which plots two states on the same grid that may be components of a coupled
        channels wave function.

        Args:
            Psi_coefficients1 (ndarray): Array of type 'double' or 'complex'
                coefficients for the representation of first state :math:`\Psi_1`
                on the FEM-DVR grid
            Psi_coefficients2 (ndarray): Array of type 'double' or 'complex'
                coefficients for the representation of first state :math:`\Psi_2`
                on the FEM-DVR grid
            plot_title_string (string): Title of plot, defaults to "Plot of
                FEM-DVR representation"
            N_plot_points (int): Number of points to plot, default 500
            make_plot (bool): Boolean that turns off/on this plotting feature,
                default true

        Returns:
            x_Plot, Psi_plot (ndarray, ndarrya): x and y coordinates of the graph
        """
        x_grid = self.x_pts
        w_grid = self.w_pts
        FEM_boundaries = self.FEM_boundaries
        N_order = self.n_order
        N_elements = len(FEM_boundaries) - 1
        dx = (FEM_boundaries[N_elements] - FEM_boundaries[0]) / (N_plot_points - 1)
        Psi1_Plot = []
        Psi2_Plot = []
        x_Plot = []
        for j in range(0, N_plot_points):
            xval = j * dx + FEM_boundaries[0]
            # rounding error can put the last point beyond end of grid, causing problems
            # for higher order Lagrange interpolation in high order DVRs
            if xval > FEM_boundaries[N_elements]:
                xval = FEM_boundaries[N_elements]
            x_Plot.append(xval)
            Psi1_Plot.append(self.psi_evaluate(xval, Psi_coefficients1))
            Psi2_Plot.append(self.psi_evaluate(xval, Psi_coefficients2))
        if make_plot:
            figure = plt.figure()
            plt.suptitle(plot_title_string, fontsize=14, fontweight="bold")
            string11 = "Re(psi_1(t))"
            string12 = "Im(psi_1(t))"
            string13 = "Abs(psi_1(t))"
            string21 = "Re(psi_2(t))"
            string22 = "Im(psi_2(t))"
            string23 = "Abs(psi_2(t))"
            plt.plot(x_Plot, np.real(Psi1_Plot), "-r", label=string11)
            plt.plot(x_Plot, np.imag(Psi1_Plot), "-g", label=string12)
            plt.plot(x_Plot, np.real(Psi2_Plot), "-r", linestyle="--", label=string21)
            plt.plot(x_Plot, np.imag(Psi2_Plot), "-g", linestyle="--", label=string22)
            plt.plot(x_Plot, np.abs(Psi1_Plot), "-k", label=string13)
            plt.plot(x_Plot, np.abs(Psi2_Plot), "-k", linestyle="--", label=string23)
            plt.legend(loc="best")
            plt.xlabel(" x ", fontsize=14)
            plt.ylabel("psi", fontsize=14)
            #     limits if necessary
            # xmax = float(rmax)  # CWM: need to use float() to get plt.xlim to work to set x limits
            # plt.xlim([0,xmax])
            # rmax=16.0
            # x_max_plot = float(rmax)  # CWM: need to use float() to get plt.xlim to work to set x limits
            # plt.xlim([0,x_max_plot])
            print(
                "\n Running from terminal, close figure window to proceed and make .pdf file of figure"
            )
            plt.savefig("Plot_Output/" + plot_title_string + ".pdf", transparent=False)
            plt.show()  # note plt.show() evidently clears everything for this plot
        return (
            x_Plot,
            Psi1_Plot,
            Psi2_Plot,
        )  # returns the x, y coordinates for a graph

    def Crank_Nicolson(
        self,
        t_initial,
        t_final,
        N_times,
        Coefs_at_t_initial,
        potential,
        Is_H_time_dependent=True,
    ):
        """
        Crank Nicolson Propagator for time-dependent or time-independent Hamiltonian
        Both implementations are unitary.  A time step is

        .. math::

            C_t = \Big(1 + i H\\big(t-\\frac{\Delta t}{2}\\big) *\\frac{\Delta
            t}{2}\Big)^{-1} * \Big(1 - i H\\big(t-\\frac{\Delta t}{2}\\big)
            *\\frac{\Delta t}{2}\Big)*C_t-\Delta t

        The ECS FEM-DVR Hamiltonian is complex symmetric, and is not conjugated in
        this expression

        Args:
            t_initial (int): Initial time for Crank_Nicolson propagation
            t_final (int): Final time for Crank_Nicolson propagation
            N_times (int): Number of time steps
            Coefs_at_t_initial (ndarray): Array of coefficients representing the
                wavefunction in the FEM-DVR grid at initial time
            potential (function): A caller provided potential that provides the
                pertubation to the system
            Is_H_time_dependent (bool): Boolean for specifying if Hamiltonian is
            time-dependent, defaults to true

        Returns:
            Ct (ndarray): new array of coefficients after one Crank Nicolson time-step
        """
        nbas = self.nbas
        x = self.x_pts
        t_interval = t_final - t_initial
        Deltat = t_interval / np.float(N_times)
        print(
            "Start Crank Nicolson propagation from ",
            t_initial,
            " to ",
            t_final,
            " atu" "\n                              in steps of ",
            Deltat,
            " with ",
            N_times,
            " steps",
        )
        # Extract the kinetic energy matrix calculated in earlier call to DVRHelper()
        # KE = np.zeros((nbas,nbas), dtype=np.complex)
        KE = np.copy(self.KE_mat)
        # Copy Coefs_at_t_initial, so they won't be changed upon return
        CPrevious = np.zeros(nbas, dtype=np.complex)
        for i in range(0, nbas):
            CPrevious[i] = Coefs_at_t_initial[i]
        # initial construction of the matrices of (1-i H Deltat/2) and (1+i H Deltat/2)
        # both with H(t = t_initial)
        Ct = np.zeros(nbas)
        M = np.identity(nbas, dtype=np.complex)
        M = M + 1j * KE * Deltat / 2.0
        Mconj = np.identity(nbas, dtype=np.complex)
        Mconj = Mconj - 1j * KE * Deltat / 2.0
        # vectorized logic for potential "correction" meaning diagonal part of M matrices
        i = np.arange(nbas)
        potential_correction = (1j * Deltat / 2.0) * potential(x[i + 1], t_initial)
        # ordinary for loop (slower) for potential correction
        # to use instead if vectorization of potential() gives errors at execution like:
        #  in _vectorize_call res = array(outputs, copy=False, subok=True, dtype=otypes[0])
        #  TypeError: can't convert complex to float
        #
        #        potential_correction  = np.zeros(nbas, dtype=np.complex)
        #        for i in range(0,nbas):
        #            potential_correction[i] =  (1j*Deltat/2.0)*potential(x[i+1],t_initial)
        #
        np.fill_diagonal(M, M.diagonal() + potential_correction)
        np.fill_diagonal(Mconj, Mconj.diagonal() - potential_correction)
        #
        if Is_H_time_dependent:
            #  open for loop on time steps, starting at t = Deltat
            for itime in range(1, N_times + 1):
                time = itime * Deltat + t_initial
                #  before each step update the diagonals (only, for efficiency) of M and M* matrices

                potential_correction_M = 1.0 + (1j * Deltat / 2.0) * (
                    KE[i, i] + potential(x[i + 1], time - Deltat / 2.0)
                )
                potential_correction_Mconj = 1.0 - (1j * Deltat / 2.0) * (
                    KE[i, i] + potential(x[i + 1], time - Deltat / 2.0)
                )
                np.fill_diagonal(M, potential_correction_M)
                np.fill_diagonal(Mconj, potential_correction_Mconj)
                #  step to time t from t - Delta t is
                #  C_t = (1 + i H(t-Deltat/2)*Deltat/2)^-1 *(1 - i H(t-Deltat/2)*Deltat/2) C_t-Deltat
                #  updated  M before making M*
                rhs = Mconj.dot(CPrevious)
                Ct = LA.solve(M, rhs)
                CPrevious = np.copy(Ct)
        else:
            # compute propagation matrix U = M^-1 M*
            Minv = LA.inv(M)
            U = np.matmul(Minv, Mconj)
            #  open for loop on time steps, starting at t = Deltat
            for itime in range(1, N_times + 1):
                time = itime * Deltat + t_initial
                #  step to time t from t - Delta t is
                #  C_t = (1 + i H(t)*Deltat/2)^-1 *(1 - i H(t-Deltat)*Deltat/2) C_t-Deltat
                #  so for time-independent H
                #  C_t = U*C_t-Deltat
                Ct = U.dot(CPrevious)
                CPrevious = np.copy(Ct)

        print("End Crank-Nicolson propagation at time = ", time)
        return Ct

    def Crank_Nicolson_Two_States(
        self,
        t_initial,
        t_final,
        N_times,
        Coefs_at_t_initial,
        V_potential_1,
        V_potential_2,
        V_coupling,
        Is_H_time_dependent=True,
    ):
        """
        Crank Nicolson Propagator for time-dependent or time-independent Hamiltonian
        Both implementations are unitary.  A time step is

        .. math::

            C_t = \Big(1 + i H\\big(t-\\frac{\Delta t}{2}\\big) *\\frac{\Delta
            t}{2}\Big)^{-1} * \Big(1 - i H\\big(t-\\frac{\Delta t}{2}\\big)
            *\\frac{\Delta t}{2}\Big)*C_t-\Delta t

        In this routine the wave function has two components propagating on
        coupled potential curves specified by V_potential_1, V_potential_2 and V_coupling

        Time-dependent case can slow.  Building the Hamiltonian and M matrices
        dominates for small grids, but linear equation solve dominates time for
        large grids. This implementation allows both the diagonal and coupling
        potentials to be time-dependent and allows the coupling to be complex hermitian.

        Args:
            t_initial (int): Initial time for Crank_Nicolson propagation
            t_final (int): Final time for Crank_Nicolson propagation
            N_times (int): Number of time steps
            Coefs_at_t_initial (ndarray): Array of coefficients representing the
                wavefunction in the FEM-DVR grid at initial time
            V_potential_1 (function): A caller provided potential for the
                propagation of the first component of the wavefunction
            V_potential_2 (function): A caller provided potential for the
                propagation of the second component of the wavefunction
            V_coupling (function): A caller provided potential that provides
                coupling between the two potential curves
            Is_H_time_dependent (bool): Boolean for specifying if Hamiltonian is
                time-dependent, defaults to true

        Returns:
            Ct (ndarray): new array of coefficients after one Crank Nicolson time-step
        """
        nbas = self.nbas
        t_interval = t_final - t_initial
        Deltat = t_interval / np.float(N_times)
        print(
            "Start Crank Nicolson propagation from ",
            t_initial,
            " to ",
            t_final,
            " atu" "\n                              in steps of ",
            Deltat,
            " with ",
            N_times,
            " steps",
        )
        # Build Hamiltonian at t_initial
        H_mat = self.Hamiltonian_Two_States(
            V_potential_1, V_potential_2, V_coupling, t_initial
        )
        # Copy Coefs_at_t_initial, so they won't be changed upon return
        CPrevious = np.zeros(2 * nbas, dtype=np.complex)
        for i in range(0, 2 * nbas):
            CPrevious[i] = Coefs_at_t_initial[i]
        # initial construction of the matrices of (1-i H Deltat/2) and (1+i H Deltat/2)
        # both with H(t = t_initial)
        Ct = np.zeros(2 * nbas)
        M = 1j * H_mat * Deltat / 2.0
        # not built with complex conj. because H_mat can be complex hermition
        Mconj = -1j * H_mat * Deltat / 2.0
        potential_correction = 1.0
        np.fill_diagonal(M, M.diagonal() + potential_correction)
        np.fill_diagonal(Mconj, Mconj.diagonal() + potential_correction)
        #
        #  Start the propagation, checking flag for time-dependent Hamiltonian
        if Is_H_time_dependent:
            #  open for loop on time steps, starting at t = Deltat
            for itime in range(1, N_times + 1):
                clock_start = timeclock.time()
                time = itime * Deltat + t_initial
                # build full H at time-Deltat/2, for step from t-Deltat to t
                H_mat = self.Hamiltonian_Two_States(
                    V_potential_1, V_potential_2, V_coupling, time - Deltat / 2.0,
                )
                Ct = np.zeros(2 * nbas)
                M = 1j * H_mat * Deltat / 2.0
                Mconj = -1j * H_mat * Deltat / 2.0
                potential_correction = 1.0
                np.fill_diagonal(M, M.diagonal() + potential_correction)
                np.fill_diagonal(Mconj, Mconj.diagonal() + potential_correction)
                rhs = Mconj.dot(CPrevious)
                #  clock_build = timeclock.time() # uncomment timing logic for detailed timings
                Ct = LA.solve(M, rhs)
                CPrevious = np.copy(Ct)
                #  clock_solve = timeclock.time()
                #  print("time for H and M build = ",clock_build - clock_start, \
                #            " time for solve = ",clock_solve-clock_build," secs")
        else:
            # compute propagation matrix U = M^-1 M* which is used for all subsequent steps
            Minv = LA.inv(M)
            U = np.matmul(Minv, Mconj)
            #  open for loop on time steps, starting at t = Deltat
            for itime in range(1, N_times + 1):
                time = itime * Deltat + t_initial
                #  step to time t from t - Delta t is
                #  C_t = (1 + i H(t)*Deltat/2)^-1 *(1 - i H(t-Deltat)*Deltat/2) C_t-Deltat
                #  so for time-independent H
                #  C_t = U*C_t-Deltat
                Ct = U.dot(CPrevious)
                CPrevious = np.copy(Ct)

        print("End Crank-Nicolson propagation at time = ", time)
        return Ct

    def rescale_kinetic(self, reduced_mass):
        """
        initialize KE_mat with the Kinetic energy matrix multiplied by
        :math:`\\frac{1}{\mu}`, where :math:`\mu` is the reduced mass

        Args:
            reduced_mass (double): The reduced mass used to rescale the KE
                matrix

        Returns:
            KE_mat (ndarray): rescaled KE matrix
        """
        self.KE_mat = (1.0 / reduced_mass) * self.KE_mat

    def hello_world(self):
        print("HelloWorld from FEM_DVR class!")
