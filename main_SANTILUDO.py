"""
PROGRAM SANTILUDO
- Programmers: *S.G. Solazzi & L. Bodet* (vanGen.m courtesy of D. Jougnot)
- Frist version: *2020/07/14*
- Last update (first draft with synthetic computation) on 2022/03
- Some minor modifs 2024/03
- Converted to C++, and interfaced with Python, on 2024/03 by J. Cunha Teixeira
- Adapted to multi-layers, on 2024/04 by J. Cunha Teixeira

## Objectives:
This program takes the hydromechanical properties of a partially 
saturated 1D soil for a given steady-state condition and computes:
- (i)   The saturation profile in depth ;
- (ii)  The variations of VPs, VSs, and rhobs with saturation (depth) ;    
- (iii) The P-wave first arrival times and suface-wave dispersion 
         for a given linear acquisition setup.

## Reference
Solazzi, S. G., Bodet, L., Holliger, K., & Jougnot, D. (2021). Surface-wave 
dispersion in partially saturated soils: The role of capillary forces. Journal
of Geophysical Research: Solid Earth, 126, e2021JB022074. 
https://doi.org/10.1029/2021JB022074


## Dependencies:
  Source code is stored in ./src folder.
  Functions are written in C++ file (*.cpp) and wrapped in Cython file (*.pyx) to be compiled
  as a shared library and called in Python.
  - Source functions in C++ are named as nameFunction_src
  - Wapped functions in Cython to be called in Python are named as nameFunction
  -> If you change the header of a function in the *.cpp, you must change it on the *.pyx file also.

  ### C++ functions for RP and VG models are stored in ./src/RPfunctions_src.cpp and ./src/VGfunctions_src.cpp
      - ./src/GVfunctions_src.cpp/vanGen_src -> to compute Saturation from...
      - ./src/GVfunctions_src.cpp/hillsAverage_src -> to compute ... from ...
      - ...
  ### Seismic data forward modelling utils stored in ./src/TTDSPfunctions_src.cpp :
      - ./src/TTDSPfunctions_src.cpp/firstArrival_src -> compute P- or S-wave first arrival times from a 1D velocity model
      + *gpdc* -> code (to be compiled and called) from https://www.geopsy.org to to compute surface-wave dispersion from a 1D velocity model
      + I/O functions:
      - ./src/TTDSPfunctions_src.cpp/writeVelocityModel_src -> create 1D velocity models readable by *gpdc*
      - ./src/TTDSPfunctions_src.cpp/readDispersion -> Only function not wtitten in C++ but directly in Cython in the *.pyx file
        reads dispersion curves created by *gpdc*
 ### Utils in bash and or latex to plot results ./src folder:
      - *cropSANTILUDO.sh* (as to be chmod +x...) 

## Compilation is done with the setup.py file :
python3 setup.py build_ext --inplace clean
This creates a ./build and ./bib directories, with the last containing the shared libraries *.so
that can be used in Python and the resulting *.cpp files from the *.pyx compilation
"""





# Python packages
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import copper
from time import perf_counter
from io import StringIO
from subprocess import run, PIPE

# Cython shared libraries for C++ and Cython functions
from lib.RPfunctions import *
from lib.VGfunctions import *
from lib.TTDSPfunctions import *



# Ploting font size parameter definition
plt.rcParams.update({'font.size': 18})


# Create output folder if do not exist
if not os.path.exists("./output/"):
   os.makedirs("./output/")



### TIME START --------------------------------------------------------------------------------------------------------------------------------------
start = perf_counter()
### -------------------------------------------------------------------------------------------------------------------------------------------------



### ROCK PHYSICS PARAMETERS -------------------------------------------------------------------------------------------------------------------------
# General physical constants
rhow = 1000.0 # Water density [Kg/m3]
rhoa = 1.0 # Air density [Kg/m3]
kw = 2.3e9 # Water bulk modulus [Pa]
ka = 1.01e5 # Air bulk modulus [Pa]
g = 9.82 # Gravity acceleration [m/s2]

# Soil layers
soiltypes = ['sand', 'clay',] # Soil types (see list in selectSoilType.m with associated VG parameters and mixture)
thicknesses = [5, 5] # Layer thicknesses [m]

# Grains/agregate parameters per layer
Ns = [8, 8] # Coordination Number (number of contact per grain) | default = 8
fracs = [0.3, 0.3] # Fraction of non-slipping grains (helps making the soil less stiff)

# Geometry and discretisation of the medium
depth = np.sum(thicknesses) # Depth of the soil column [m]
dz = 0.01 # Depth sample interval [m]
top_surface_level = dz # Altitude of the soil surface[m]
zs = -np.arange(top_surface_level, depth + dz, dz) # Depth positions (negative downward) [m]
NbCells = len(zs) - 1 # Number of exploration points in depth [#]

# Water table
WTs = [2, 7] # Water table depths [m]
color_map = copper(np.linspace(0, 1, len(WTs))) # Colorscale for plots if several WT tested

# Grains/agregate mechanical properties
mu_clay = 6.8 # Shear moduli [GPa]
mu_silt = 45.0
mu_sand = 45.0
k_clay = 25.0 # Bulk moduli [GPa]
k_silt = 37.0
k_sand = 37.0
rho_clay = 2580.0 # Density [kg/m3]
rho_silt = 2600.0
rho_sand = 2600.0

# Three possible RP models:
# kk = 1 # Constant Pe (see the approach of Zyserman et al., 2017)
# kk = 2 # Pe without suction
kk = 3 # Pe with suction (cf. Solazzi et al. 2021)
### -------------------------------------------------------------------------------------------------------------------------------------------------



### SEISMIC PARAMETERS ------------------------------------------------------------------------------------------------------------------------------
# Layers to put under the studied soil column on the velocity model
# In GPDC format : "thickness Vp Vs rho\n"
# Each layer is separated by \n | Only spaces between values | Last layer thickness must be 0)
# under_layers = "" # Empty string if no under layers
under_layers = "0 2400 1200 2500\n" # One substratum layer
n_under_layers = under_layers.count('\n') # Number of under layers


thks = np.diff(np.abs(zs)) # thickness vector [m]


x0 = 1 # first geophone position [m]
Nx = 192 # number of geophones [m]
dx = 1 # geophone interval [m]
xs = np.arange(x0, Nx * dx + 1, dx)
trig  = 0 # data pretrig (if needed)


# Frequency domain and sampling setup to compute dispersion
nf = 46 # number of frequency samples [#]
df = 1 # frequency sample interval [Hz]
min_f = 5 # minimum frequency [Hz]
max_f = min_f + (nf - 1) * df

n_modes = 1 # Number of modes to compute
s = 'frequency' # Over frequencies mode
wave = 'R' # Rayleigh (PSV) fundamental mode
### -------------------------------------------------------------------------------------------------------------------------------------------------



### CHECK PARAMETERS --------------------------------------------------------------------------------------------------------------------------------
if len(set(map(len, (soiltypes, thicknesses, Ns, fracs)))) != 1:
    raise ValueError(f"Arrays are not the same size : {soiltypes = }, {thicknesses = }, {Ns = }, {fracs = }")


layers = [layer for layer in under_layers.split('\n') if layer.strip()]
for i, layer in enumerate(layers) :
    if len(np.fromstring(layer, dtype=float, sep=' ')) != 4 :
        raise ValueError(f"Under layers format is not correct:\nLayer {i+1} is '{layer}'\nCorrect format : 'thickness Vp Vs rho'")
    if i == len(layers) - 1 and np.fromstring(layer, dtype=float, sep=' ')[0] != 0 :
        raise ValueError(f"Last layer thickness must be 0\nLayer {i+1} is '{layer}'\nCorrect format : '0 Vp Vs rho'")
    if i < len(layers) - 1 and np.fromstring(layer, dtype=float, sep=' ')[0] == 0 :
        raise ValueError(f"Layers thickness must be greater than 0 (except for the last layer)\nLayer {i+1} is '{layer}'\nCorrect format : 'thickness Vp Vs rho'")
### -------------------------------------------------------------------------------------------------------------------------------------------------



# Runnung the program for each water table depth
for iWT, WT in enumerate(WTs) :
    #### ROCK PHYSICS -------------------------------------------------------------------------------------------------------------------------------
    # Saturation profile with depth
    hs, Sws, Swes = vanGen(zs, WT, soiltypes, thicknesses)
    # fig, axs = plt.subplots(2, 1)
    # fig.suptitle("vanGen")
    # axs[0].plot(hs, zs)
    # axs[0].axhline(-WT, color='gray', linestyle='--')
    # axs[0].set_xlabel('Pressure head h')
    # axs[0].set_ylabel('Depth z [m]')
    # axs[1].plot(Sws, zs, Swes, zs)
    # axs[1].axhline(-WT, color='gray', linestyle='--')
    # axs[1].set_xlabel('Saturation S')
    # axs[1].set_ylabel('Depth z [m]')
    # axs[1].legend(['Total saturation Sw', 'Effective wetting phase saturation Swe'])
    # axs[1].set_xlim([0,1.05])

    # Effective Grain Properties (constant with depth)
    mus, ks, rhos, nus = hillsAverage(mu_clay, mu_silt, mu_sand, rho_clay,
                                      rho_silt, rho_sand, k_clay, k_silt,
                                      k_sand, soiltypes)
    # print('\nhillsAverage')
    # print(f'Shear moduli of grains: {mus = } [Pa]')
    # print(f'Bulk moduli of grains: {ks = } [Pa]')
    # print(f'Density of  grains: {rhos = } [kg/m3]')
    # print(f"Poisson's ratio: {nus = }")
    # print('\n')

    # Effective Fluid Properties
    kfs, rhofs, rhobs = effFluid(Sws, kw, ka, rhow,
                                 rhoa, rhos, soiltypes, thicknesses, dz)
    # fig, axs = plt.subplots(2, 1)
    # fig.suptitle("effFluid")
    # axs[0].plot(kfs, zs)
    # axs[0].axhline(-WT, color='gray', linestyle='--')
    # axs[0].set_xlabel('Effective compressibility k_f [Pa-1]')
    # axs[0].set_ylabel('Depth z [m]')
    # axs[1].plot(rhofs, zs, rhobs, zs)
    # axs[1].axhline(-WT, color='gray', linestyle='--')
    # axs[1].set_xlabel('Density rho [kg/m3]')
    # axs[1].set_ylabel('Depth z [m]')
    # axs[1].legend(['Effective fluid density rhof', 'Bulk density rhob'])

    # Hertz Mindlin Frame Properties
    KHMs, muHMs = hertzMindlin(Swes, zs, hs, rhobs,
                               g, rhoa, rhow, Ns,
                               mus, nus, fracs, kk,
                               soiltypes, thicknesses)
    # fig, ax = plt.subplots()
    # fig.suptitle("hertzMindlin")
    # ax.plot(KHMs, zs)
    # ax.plot(muHMs, zs)
    # ax.axhline(-WT, color='gray', linestyle='--')
    # ax.set_xlabel('Pressure [Pa]')
    # ax.set_ylabel('Depth z [m]')
    # ax.legend(['Effective bulk KHM', 'Shear moduli muHM'])
    # ax.set_xlim([1.3e8,2.5e8])

    # Saturated Properties
    VPs, VSs = biotGassmann(KHMs, muHMs, ks, kfs,
                            rhobs, soiltypes, thicknesses, dz)
    # fig, ax = plt.subplots()
    # fig.suptitle("biotGassmann")
    # ax.plot(VPs, zs)
    # ax.plot(VSs, zs)
    # ax.axhline(-WT, color='gray', linestyle='--')
    # ax.set_xlabel('Velocity [m/s]')
    # ax.set_ylabel('Depth z [m]')
    # ax.legend(['Vp', 'Vs'])

    # plt.show()
    # plt.close()
    ### ---------------------------------------------------------------------------------------------------------------------------------------------



    #### SEISMIC FWD MODELING -----------------------------------------------------------------------------------------------------------------------
    # First arrival time computations
    ThodPs = firstArrival(thks, VPs, xs, trig) # P-wave first arrival times
    ThodSs = firstArrival(thks, VSs, xs, trig) # S-wave first arrival times

    # Velocity model in string format for GPDC
    velocity_model_string = writeVelocityModel(thks, VPs, VSs, rhobs, under_layers, n_under_layers)

    # Dispersion curves computing with GPDC
    velocity_model_RAMfile = StringIO(velocity_model_string) # Keep velocity model string in the RAM in a file format alike to trick GPDC which expects a file
    gpdc_command = [f"gpdc -{wave} {n_modes} -n {nf} -min {min_f} -max {max_f} -s {s}"]
    gpdc_output_string = run(gpdc_command, input=velocity_model_RAMfile.getvalue(), text=True, shell=True, stdout=PIPE).stdout # raw output string from GPDC

    dispersion_data, n_modes = readDispersion(gpdc_output_string) # Reads GPDC output and converts dispersion data to a list of numpy arrays for each mode
                                                                  # Updates number of computed modes (can be lower than what was defined if frequency range too small)
    ### ---------------------------------------------------------------------------------------------------------------------------------------------



    ### PLOTS ---------------------------------------------------------------------------------------------------------------------------------------
    color = color_map[iWT,:] # Color for the plot

    if iWT == 0:
        # Create figures
        fig1, ax1 = plt.subplots(dpi=300)
        fig2, ax2 = plt.subplots(dpi=300)
        fig3, ax3 = plt.subplots(dpi=300)
        fig4, ax4 = plt.subplots(dpi=300)
        fig5, ax5 = plt.subplots(dpi=300)
        fig1000, ax1000 = plt.subplots(dpi=300)
        fig2000, ax2000 = plt.subplots(dpi=300)
        fig3000, ax3000 = plt.subplots(dpi=300)

        # Name for the output files with soil types and thicknesses
        name = ''
        for soiltype, thickness in zip(soiltypes, thicknesses):
            name += f"{soiltype}{thickness}_"
        name = name[:-1]

    # Plot Sws vs. zs
    ax1.plot(Sws, zs, linewidth=2, color=color)
    ax1.axhline(-WT, color=color, linestyle='--', linewidth=1)
    ax1.set_xlim(0, 1.1)
    ax1.set_xlabel('$S_w$')
    ax1.set_ylabel('$z$ [m]')
    fig1.savefig('./output/1.pdf', bbox_inches='tight')

    # Plot VPs vs. zs
    ax2.plot(VPs, zs, linewidth=2, color=color)
    ax2.axhline(-WT, color=color, linestyle='--', linewidth=1)
    ax2.set_xlim(0, 2000)
    ax2.set_xlabel('$V_p$ [m/s]')
    ax2.set_ylabel('$z$ [m]')
    fig2.savefig('./output/2.pdf', bbox_inches='tight')

    # Plot VSs vs. zs
    ax3.plot(VSs, zs, linewidth=2, color=color)
    ax3.axhline(-WT, color=color, linestyle='--', linewidth=1)
    ax3.set_xlim(0, 600)
    ax3.set_xlabel('$V_s$ [m/s]')
    ax3.set_ylabel('$z$ [m]')
    fig3.savefig('./output/3.pdf', bbox_inches='tight')

    # Plot rho_b vs. zs
    ax4.plot(rhobs, zs, linewidth=2, color=color)
    ax4.axhline(-WT, color=color, linestyle='--', linewidth=1)
    ax4.set_xlim(1500, 2000)
    ax4.set_xlabel('$\\rho_b$ [kg/$m^3$]')
    ax4.set_ylabel('$z$ [m]')
    fig4.savefig('./output/4.pdf', bbox_inches='tight')

    # Plot Poisson ratio vs. zs
    ax5.plot(list(map(fish, VPs, VSs)), zs, linewidth=2, color=color)
    ax5.axhline(-WT, color=color, linestyle='--', linewidth=1)
    ax5.set_xlim(0, 0.5)
    ax5.set_xlabel('Poisson ratio')
    ax5.set_ylabel('$z$ [m]')
    fig5.savefig('./output/5.pdf', bbox_inches='tight')

    # Plot simulated P-wave first arrivals
    ax1000.plot(xs, ThodPs, linewidth=2, color=color)
    ax1000.set_xlim([0, max(xs)])
    ax1000.set_ylim([0, 0.25])
    ax1000.set_xlabel('Offset [m]')
    ax1000.set_ylabel('P- first arrival time [s]')
    fig1000.savefig('./output/1000.pdf', bbox_inches='tight')

    # Plot simulated S-wave first arrivals
    ax2000.plot(xs, ThodSs, linewidth=2, color=color)
    ax2000.set_xlim([0, max(xs)])
    ax2000.set_ylim([0, 0.7])
    ax2000.set_xlabel('Offset [m]')
    ax2000.set_ylabel('S- first arrival time [s]')
    fig2000.savefig('./output/2000.pdf', bbox_inches='tight')

    # Plot simulated dispersion curve
    max_vr = 0
    min_vr = 1e10
    for mode in range(n_modes):
        ax3000.plot(dispersion_data[mode][:,0], dispersion_data[mode][:,1], linewidth=2, color=color)
        if np.max(dispersion_data[mode][:,1]) > max_vr:
            max_vr = np.max(dispersion_data[mode][:,1])
        if np.min(dispersion_data[mode][:,1]) < min_vr:
            min_vr = np.min(dispersion_data[mode][:,1])
        # Save dispersion data to file
        # np.savetxt(f'./output/{name}_M{mode}_WT{WT}.txt', dispersion_data[mode], fmt='%f')
    ax3000.set_xlim([min_f-5, max_f+5])
    ax3000.set_ylim([min_vr-100, max_vr+100])
    ax3000.set_xlabel('Frequency [Hz]')
    ax3000.set_ylabel('P-SV phase vel. [m/s]')
    fig3000.savefig('./output/3000.pdf', bbox_inches='tight')
    ### ---------------------------------------------------------------------------------------------------------------------------------------------



### CONVERT PLOTS TO SINGLE PDF ---------------------------------------------------------------------------------------------------------------------
plt.close('all')
# Change directory
os.chdir('./src')
# Change file permissions
run(['chmod', '+x', 'cropSANTILUDO.sh'])
# Execute script with arguments
run(['sh', 'cropSANTILUDO.sh', name])
# Change back to previous directory
os.chdir('..')
### -------------------------------------------------------------------------------------------------------------------------------------------------



### TIME END ----------------------------------------------------------------------------------------------------------------------------------------
end = perf_counter()
print(f"\nElapsed time : {end-start} s\n")
### -------------------------------------------------------------------------------------------------------------------------------------------------