"""
Class for a more OOP interface to several potentials.
"""
# Import  NumPy which is used to define pi, sqrt, array, .transpose etc. as
import numpy as np

from scipy.interpolate import CubicSpline


class Potential(object):
    def __init__(self, file=None):
        """Constructor method, added vectorized version of all the methods.

        Args:
            file (string, optional): The caller may provide a file to interpolate a
             potential onto the dvr grid. Default to None

        Attributes:
            vectorized_V_morse (ndarray): Vectorized version of the morse
             function
            vectorized_V_Bernstein (ndarray): Vectorized version of the
            Bernstein function
            vectorized_V_c_state (ndarray): Vectorized version of the c-state
            interpolated function if the cStateDCalc.csv file is provided

        Note:
            At this point constants are for Lucchese 4/3/2020 calculation:
            :math:`c ^4\Sigma_u^-` state of :math:`O_2^+` where the orbitals come from a SA-MCSCF
            on the ion using an aug-cc-vTZP basis set.

        Todo:
            Add a general interpolation scheme so any file passed into this
            class's constructor will work

        """

        self.vectorized_V_morse = np.vectorize(self.V_morse)
        self.vectorized_V_Bernstein = np.vectorize(self.V_Bernstein)
        self.vectorized_V_c_state = np.vectorize(self.V_c_state)

        if file is None:
            print(
                "Potential constructed without a file. Can only use analytic potential functions defined in the quantumgrid.potential class"
            )
        else:
            file_name_c_state = open(file, "r")
            data = np.loadtxt(file_name_c_state, delimiter=",")
            pot_len_c_state = data.shape[0]
            pot_columns = data.shape[1]
            print(
                "Finished reading file with c state potential with ",
                pot_len_c_state,
                " rows and ",
                pot_columns,
                " columns",
            )
            self.r_vals_c_state = np.empty(pot_len_c_state)
            self.V_vals_c_state = np.empty(pot_len_c_state)
            for i in range(0, pot_len_c_state):
                self.r_vals_c_state[i] = data[i, 0]
                self.V_vals_c_state[i] = data[i, 1]

    def V_morse(self, r, t=0.0):
        """
        Morse Potential defined by

        .. math::

            V = d*(y^2 - 2*y) + \\mathrm{Centrifugal potential}

        with :math:`y` defined by

        .. math::

         y = e^{(-a*(r-re))}

        This potential also defines parameters specifically for :math:`H_2`

        Args:
            r (complex): FEM-DVR point where this potential is evaluated at
            t (int): Time dependence of this potential to simulate
                        turning on a field pertubation, for example. Defaults to
                        t=0.0

        Returns:
            pot (complex): potential value at the point r at the time t
        """
        d = 0.1746
        a = 1.0277
        re = 1.4022
        H_Mass = 1.0078
        Daltons_to_eMass = 1822.89
        mu = (H_Mass / 2.0) * Daltons_to_eMass
        y = np.exp(-a * (r - re))
        # j value for centrifugal potential.  mu defined in main part of script above
        j = 0  # Morse potential has rotational predissociation resonances for some j
        pot = d * (y ** 2 - 2.0 * y) + np.float(j * (j + 1)) / (
            2.0 * mu * r ** 2
        )
        return pot

    def V_Bernstein(self, r, t=0.0):
        """
        :math:`H_2` potential from T-G. Wiechand R.B. Bernstein, J. Chem. Phys. 46 (1967) 4905.
        This is an accurate fit to the accurate Kolos and Wolneiwicz potential curve
        representation is valid :math:`0.4 <= R` to infinity
        used in old ECS calculation in
        Julia Turner and C. William McCurdy, Chemical Physics 71(1982) 127-133
        for resonances in dissociation for :math:`j .ne. 0`

        Note:
            ECS contour must begin beyond :math:`r = 9.5 a_0` for safe analytic continuation

        Args:
            r (complex): FEM-DVR point where this potential is evaluated at
            t (int): Time dependence of this potential to simulate
                        turning on a field pertubation, for example. Defaults to
                        t=0.0

        Returns:
            pot (complex): potential value at the point r at the time t
        """
        a_vec = [
            -3.7623236364e-3,
            1.4291725467e-2,
            -2.6491493104e-2,
            3.0802158643e-2,
            -2.4414431427e-2,
            1.2072690633e-2,
            1.0669803453e-2,
            -3.1351262502e-2,
            -2.4593504473e-2,
            9.0968827782e-2,
            8.0055110345e-2,
            -2.2685375608e-1,
            -1.4912492825e-1,
            3.9041633873e-1,
            1.7916153661e-1,
            -4.7291514961e-1,
            -1.4317771747e-1,
            4.1382169150e-1,
            7.3590396723e-2,
            -2.6524118029e-1,
            -1.9970631183e-2,
            1.2463802250e-1,
            -1.2491070013e-3,
            -4.2434523716e-2,
            3.4575120517e-3,
            1.0180959606e-2,
            -1.4411614262e-3,
            -1.6314090918e-3,
            3.1362830316e-4,
            1.5666712172e-4,
            -3.6848921690e-5,
            -6.8198927741e-6,
            1.8540052417e-6,
        ]
        # from  Hirshfelder and Lowdin
        # Hirshfelder and Lowdin corrected values in 1965 -1 -C6/r^6 -C8/r^8
        C6 = 6.499026
        C8 = 124.395
        # Chan and Dalgarno give their values in Rydbergs evidently. This is from paper cited by Bernstein above
        C10 = 6571.0 / 2.0
        # print("length of a_vec = ",len(a_vec))
        #
        H_Mass = 1.0078
        Daltons_to_eMass = 1822.89
        mu = (H_Mass / 2.0) * Daltons_to_eMass
        if np.real(r) >= 0.4 and np.real(r) <= 9.5:
            vsum = 0.0
            for n in range(0, 33):
                vsum = vsum + a_vec[n] * ((r - 5.0) / 2.5) ** n
                # print("n, a_vec,",n," ",'{:.10e}'.format(a_vec[n]))
        else:
            vsum = -C6 / r ** 6 - C8 / r ** 8 - C10 / r ** 10
        # j = 17 is Fig 2 of Turner+McCurdy, E_res = (0.004044878419994 -0.000219496448j)  hartrees
        j = 17
        vpot = vsum + float(j * (j + 1)) / (2.0 * mu * r ** 2)
        return vpot

    def V_c_state(self, r, t=0.0):
        """
        Interpolate computed values using scipy CubicSpline
        :math:`\\frac{1}{R^4}` tail added matching value and finite diff
        derivative at :math:`R=5`
        At this point constants are for Lucchese 4/3/2020 calculation:
        c-4-Sigma-u-(-) state of O2+ where the orbitals come from a SA-MCSCF
        on the ion using an aug-cc-vTZP basis set.

        Args:
            r_vals_c_state (ndarray): position values for potential to be
                                        interpolated
            V_vals_c_state (ndarray): Potential values interpolated on FEM-DVR
                                        grid
            r (complex): FEM-DVR point where this potential is evaluated at
            t (int): Time dependence of this potential to simulate
                        turning on a field pertubation, for example. Defaults to
                        t=0.0

        Returns:
            pot (complex): potential value at the point r at the time t
        """
        Hartree_to_eV = 27.211386245988  # NIST ref

        cs = CubicSpline(self.r_vals_c_state, self.V_vals_c_state)
        n_vals = self.r_vals_c_state.shape[0]
        if self.r_vals_c_state[0] <= r and r <= self.r_vals_c_state[n_vals - 1]:
            x = np.real(r)
            pot = cs(x)
        if np.real(r) > 5.0:
            pot = 20.26002003285 - 85.94654796874 / r ** 4
        if np.real(r) < self.r_vals_c_state[0] and np.real(r) >= 1.5:
            pot = cs(r)
        if np.real(r) < 1.5:
            print("r out of range in V_c_state ", r)
            exit()
        # interpolated value
        potential = (pot - 20.26002003285) / Hartree_to_eV
        return potential
