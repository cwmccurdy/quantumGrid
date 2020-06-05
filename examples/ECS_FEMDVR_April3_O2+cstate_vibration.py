#%%
"""                  O2+ c state predissociating resonances
                          C.W. McCurdy
                          04/04/2020
    Time-independent Exterior Complex Scaling (ECS) FEM-DVR example

  Uses ECS_DVRHelper.py class library
 
 Finite Element Method - Discrete Variable Representation (FEM-DVR)
 for 1D Schroedinger equation using Gauss-Lobatto quadrature in
 each finite element Uses class DVRHelper() to construct FEM-DVR
 points, weights and Kinetic Energy

 For time-independent potential, this example implements Exterior
 Complex Scaling on the FEM-DVR contour.  The value of R0 and the
 complex scale factor e^(I*theta) are specified.  The representation
 of the pootential must be able to be evaluated on the complex part
 of the contour.

 Example: Finds all eigenvalues of complex scaled Hamiltonian and
          plots any one of them, specified by n_Plot
"""
#%%
# preliminaries to invoke SciPy linear algebra functions 
from scipy import linalg as LA
# and NumPy which is used to define pi, sqrt, array, .transpose etc. as 
import numpy as np
import matplotlib.pyplot as plt  # import matplotlib pyplot functions
from matplotlib import animation  # for animation from same class library
import os  # functions to manipulate files and directories
from scipy.interpolate import CubicSpline  
from DVR.ECS_DVRHelper import ECS_DVRHelper  # contains Barbalinardo/McCurdy FEM-DVR and Crank Nicolson functions
import time as timeclock  # for timing parts of the calculation during debugging
#
#================ Make Directory for Plots if it's not there already =============
#  
# detect the current working directory and print it
path = os.getcwd()  
print ("The current working directory is %s" % path)  
# define the name of the directory to be created
Plot_Output = path+'/Plot_Output'
if os.path.exists(Plot_Output):
    print("Directory for wave function plots already exists",Plot_Output)
else:
    print("Attempting to create directory for wave function plots ",Plot_Output)
    try:  
        os.mkdir(Plot_Output)
    except OSError:  
        print ("Creation of the directory %s failed" % Plot_Output)
    else:  
        print ("Successfully created the directory %s " % Plot_Output)
#=====================================FEM_DVR===================================
#  Set up the FEM DVR grid given only the Finite Element boundaries and order
#  of Gauss Lobatto quadrature,  and compute the Kinetic Energy matrix for 
#  this FEM DVR for Mass = mass set in call (atomic units).
#
# Here is where the reduced mass is set  
#H_Mass = 1.007825032 #H atom atomic mass
O_Mass =  15.99491461957 # O 16 mass from NIST tables
O_18_Mass = 17.99915961286
Daltons_to_eMass = 1822.89
mu = (O_Mass*O_Mass/(O_Mass+O_Mass))*Daltons_to_eMass 
print("reduced mass mu = ",mu)
bohr_to_Angstrom = 0.529177
Hartree_to_eV = 27.211386245988 # NIST ref
eV_to_wavennumber = 8065.54393734921 # NIST ref on constants + conversions 
Hartree_to_wavenumber = 2.1947463136320e5    # value from NIST ref on constants + conversions
HartreeToKelvin = 315773;
atu_to_fs = 24.18884326509/1000
#  Set up the FEM-DVR grid
n_order = 25 # 25 seems to be enough for O2+ on 0.1 FEMs out to 7.5, used 35 to check.``
FEM_boundaries = [1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,
                  3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,
                  5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,
                  7.0,7.1,7.2,7.3,7.4,7.5,8.0, 9.0]
scale_factor = np.exp(1j*5.0*np.pi/180.0) # scale_factor = 1 Integer for no scaling
#  The only physical result that chamnges with scale factor is the quadrature
#  of the residue of the Green's function for "partial-width" type calculation of Gamma
R0 = 5.1 # scale just beyond last calculated point 
dvr = ECS_DVRHelper(n_order, FEM_boundaries,Mass=mu, Complex_scale=scale_factor, R0_scale = R0)
print("\nFEM-DVR basis of ", dvr.nbas, " functions")
#
#   Function to define potential at x and t (if potential is time-dependent)
#   goes here
#
#   Read comma separated values file from Robert Lucchese of O2^+ quartet c state
#   potential curve
file_name_c_state = open(  'cStateDCalc.csv','r')
data = np.loadtxt(file_name_c_state, delimiter=",")
pot_len_c_state = data.shape[0]
pot_columns = data.shape[1]
print("Finished reading file with c state potential with ",pot_len_c_state," rows and ",pot_columns," columns")
r_vals_c_state=np.empty(pot_len_c_state)
V_vals_c_state=np.empty(pot_len_c_state)
for i in range(0,pot_len_c_state):
    r_vals_c_state[i] = data[i,0]
    V_vals_c_state[i] = data[i,1]
#    print(r_vals_c_state[i]," ",V_vals_c_state[i])
#
cs = CubicSpline(r_vals_c_state,V_vals_c_state)
n_vals = r_vals_c_state.shape[0]
def V_c_state(r,time):
    #
    #  Interpolate computed values using scipy CubicSpline
    #  1/R^4 tail added matching value and finite diff derivative at R=5 
    #   At this point constants are for Lucchese 4/3/2020 calculation:
    #  c-4-Sigma-u-(-) state of O2+ where the orbitals come from a SA-MCSCF
    #   on the ion using an aug-cc-vTZP basis set.
    #
    #print(" r value in V_c_state ",r)
    #print("n_vals in V_c_state ",n_vals)
    if(r_vals_c_state[0] <= r and r <= r_vals_c_state[n_vals-1]):
       x = np.real(r)
       pot = cs(x)
    if(np.real(r) > 5.0):
        pot = 20.26002003285 - 85.94654796874/r**4
    if(np.real(r) < r_vals_c_state[0] and np.real(r) >= 1.5):
        pot = cs(r)
    if(np.real(r) < 1.5):
        print("r out of range in V_c_state ",r)
        exit()
    # interpolated value
    potval= (pot - 20.26002003285 )/Hartree_to_eV 
    return potval
# in order to interface correctly to the ECS_DVRHelper.py class that is vectorized
vectorized_V_c_state = np.vectorize(V_c_state)

#
#==================================================================================
#  Plot potential on the DVR grid points on which the wavefunction is defined
print("\n Plot potential ")
print("Test V",V_c_state(2.0345+0.0*1j ,.0))
#print("Test V",V_c_state(10.0+ 0.5*1j,.0))
time = 0.0
x_Plot = []
pot_Plot = []
for j in range(0,dvr.nbas):
    x_Plot.append(dvr.x_pts[j+1])
    potval = np.real(V_c_state(dvr.x_pts[j+1],time))
    pot_Plot.append(potval)
plt.suptitle('V(x) at DVR basis function nodes', fontsize=14, fontweight='bold')
string="V"
plt.plot(np.real(x_Plot),np.real(pot_Plot),'ro',label=string)
plt.plot(np.real(x_Plot),np.real(pot_Plot),'-b')
plt.legend(loc="best")
plt.xlabel(" x ", fontsize=14)
plt.ylabel("V", fontsize=14)
print("\n Running from terminal, close figure window to proceed and make .pdf file of figure")
#   Insert limits if necessary
#ymax = float(0.05)  # CWM: need to use float() to get plt.xlim to work to set x limits
#plt.ylim([-.18,ymax])
# save plot to .pdf file
plt.savefig('Plot_Output/' + 'Plot_potential' + '.pdf', transparent=False)
plt.show()
# make a plotting file of the potential at the FEM_DVR points (real parts)
file_potential = open('Potential.dat','w')
for j in range(0,dvr.nbas):
    print(np.real(x_Plot[j]),"  ",np.real(pot_Plot[j]), file=file_potential )
#
# =============Build Hamiltonian (at t=0 if time-dependent)=================================
#     Pass name of potential function explicitly here
time = 0.0
H_mat = dvr.Hamiltonian(vectorized_V_c_state, time)
print("\n Completed construction of Hamiltonian at t = 0")
#====================================================================================
#
# Find all the eigenvalues of the Hamiltonian so we can compare with known bound state energies
# or make a plot of the spectrum -- For a time-independent Hamiltonian example here
#
# LA.eigvalsh(H_mat) is for Hermitian matrices.  LA.eigvals is for general matrices
#EigenVals = LA.eigvals(H_mat)
print("Calculating ",dvr.nbas," eigenvalues and eigenvectors for plotting eigenfunctions")
EigenVals, EigenVecs  = LA.eig(H_mat,right=True,homogeneous_eigvals=True)
print("after LA.eig()")
#
n_energy = dvr.nbas 
file_opened = open('Spectrum_ECS.dat','w')
print("EigenVals shape ",EigenVals.shape)
for  i in range(0,n_energy):
     print("E( ",i,") =   ",EigenVals[0,i]," hartrees" )
     print(np.real(EigenVals[0,i]),"  ",np.imag(EigenVals[0,i]), file=file_opened )
#====================================================================================
#
# Extract the n_Plot'th eigenfunction for plotting 
#
# pick one of the bound states of Morse Potential to plot
# numbering can depend on numpy and python installation that determines
# behavior of the linear algebra routines.
n_Plot = n_energy -1
n_Plot =  1353 # 1351 # 1353 # 1356 last two are below barrier
wfcnPlot = []
for j in range(0,dvr.nbas):
    wfcnPlot.append(EigenVecs[j,n_Plot])
#
# normalize  wave function from diagonalization
# using integration of square on the contour
# following the old Rescigno McCurdy idea for partial widths
norm_squared = 0.
for j in range(0,dvr.nbas):
    norm_squared = norm_squared + (wfcnPlot[j])**2    
wfcnPlot = wfcnPlot/np.sqrt(norm_squared)
norm_squared = 0.
gamma_residue = 0.0
k_momentum = np.sqrt(2.0*mu*np.real(EigenVals[0,n_Plot]))  # background momentum defined with Re(Eres)
#k_momentum = np.real(np.sqrt(2.0*mu*EigenVals[0,n_Plot]))
for j in range(0,dvr.nbas):
    norm_squared = norm_squared + (wfcnPlot[j])**2    
    if(dvr.x_pts[j+1] < 5.8 ):
          free_wave = (2.0*np.sqrt(mu/k_momentum))*np.sin(k_momentum*dvr.x_pts[j+1])
          gamma_residue = gamma_residue + wfcnPlot[j]*V_c_state(dvr.x_pts[j+1],time)* \
                            free_wave*np.sqrt(dvr.w_pts[j+1])
print("Complex symmetric inner product (psi|psi) is being used")
print("Norm of wave function from int psi^2 on contour being plotted is ", np.sqrt(norm_squared))
print(" For this state the asymptotic value of k = ",k_momentum)
print("gamma from int = ",gamma_residue," |gamma|^2 = ",np.abs(gamma_residue)**2)
# Plot wave function -- It must be type np.complex
Cinitial  = np.zeros((dvr.nbas), dtype=np.complex) 
wfcnInitialPlot  = np.zeros((dvr.nbas), dtype=np.complex) 
for j in range(0,dvr.nbas):
   Cinitial[j] = wfcnPlot[j]
#
# plot n_Plot'th eigenfunction 
#
print("\n Plot Hamiltonian eigenfunction number ",n_Plot," with energy ",EigenVals[0,n_Plot])
tau = atu_to_fs/(-2.0*np.imag(EigenVals[0,n_Plot]))
print("  Lifetime tau = 1/Gamma = ",tau," fs")
number_string = str(n_Plot)
title = 'Wavefunction number = '+number_string
wfn_plot_points = 2000
#x_Plot_array, Psi_plot_array = dvr.Plot_Psi(Cinitial, plot_title_string=title,N_plot_points=wfn_plot_points,make_plot=True)
x_Plot_array, Psi_plot_array = dvr.Plot_Psi(Cinitial, plot_title_string=title,N_plot_points=wfn_plot_points,make_plot=True)
#
#  Make data file for n_Plot'th eigenfunction
#  Psi and Psi**2
#
filename='wavefunction'+number_string+'.dat'
file_opened = open(filename,'w')
print_points = len(x_Plot_array)
print("x_Plot_array shape ",print_points)
for  i in range(print_points):
     free_wave = (2.0*np.sqrt(mu/k_momentum))*np.sin(k_momentum*x_Plot_array[i])
     integrand = Psi_plot_array[i]*V_c_state(x_Plot_array[i],time)*free_wave  # for partial width gamma
     #integrand = Psi_plot_array[i]**2 # for normalization integral with c-product
     print(np.real(x_Plot_array[i]),"  ",np.imag(x_Plot_array[i]),"  ",
       np.real(Psi_plot_array[i]),"  ",np.imag(Psi_plot_array[i]),"  ",
#       np.real(Psi_plot_array[i]**2),"  ",np.imag(Psi_plot_array[i]**2),
       np.real(integrand),np.imag(integrand), file=file_opened )
#
exit()


