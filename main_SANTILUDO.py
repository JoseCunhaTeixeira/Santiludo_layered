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
from matplotlib.backends.backend_pdf import PdfPages
from time import perf_counter
from io import StringIO
from subprocess import run, PIPE, CalledProcessError

# Cython shared libraries for C++ and Cython functions
from lib.VGfunctions import vanGen
from lib.RPfunctions import hillsAverage, effFluid, hertzMindlin, biotGassmann, fish
from lib.TTDSPfunctions import firstArrival, writeVelocityModel, readDispersion



# Ploting font size parameter definition
plt.rcParams.update({'font.size': 14})


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

# Soil layers
soiltypes = ['clay'] # Soil types (see list in selectSoilType.m with associated VG parameters and mixture)
thicknesses = [5] # Layer thicknesses [m]

# Grains/agregate parameters per layer
Ns = [9] # Coordination Number (number of contact per grain) | default = 8
fracs = [0.3] # Fraction of non-slipping grains (helps making the soil less stiff) | default = 0.3

# Water table
WTs = [2] # Water table depths [m]
color_map = copper(np.linspace(0, 1, len(WTs))) # Colorscale for plots if several WT tested

# Geometry and discretisation of the medium
depth = np.sum(thicknesses) # Depth of the soil column [m]
dz = 0.01 # Depth sample interval [m]
top_surface_level = dz # Altitude of the soil surface[m]
zs = -np.arange(top_surface_level, depth + dz, dz) # Depth positions (negative downward) [m]
NbCells = len(zs) - 1 # Number of exploration points in depth [#]

# Three possible RP models:
# kk = 1 # Constant Pe (see the approach of Zyserman et al., 2017)
# kk = 2 # Pe without suction
kk = 3 # Pe with suction (cf. Solazzi et al. 2021)
### -------------------------------------------------------------------------------------------------------------------------------------------------



### SEISMIC PARAMETERS ------------------------------------------------------------------------------------------------------------------------------
# Layers to put under the studied soil column on the velocity model
# In GPDC format : [thickness Vp Vs rho]
# under_layers = [] # Empty list if no under layers
under_layers = [
                [10, 4000, 2000, 2500],
                [0, 8000, 4000, 2500],
                ]
under_layers = np.array(under_layers)
n_under_layers = len(under_layers) # Number of under layers


thks = np.diff(np.abs(zs)) # thickness vector [m]


x0 = 0.5 # first geophone position [m]
Nx = 400 # number of geophones [m]
dx = 0.5 # geophone interval [m]
xs = np.arange(x0, Nx * dx + 1, dx)
trig  = 0 # data pretrig (if needed)


# Frequency domain and sampling setup to compute dispersion
nf = 36 # number of frequency samples [#]
df = 1 # frequency sample interval [Hz]
min_f = 15 # minimum frequency [Hz]
max_f = min_f + (nf - 1) * df


n_modes = 1 # Number of modes to compute
s = 'frequency' # Over frequencies mode
wave = 'R' # Rayleigh (PSV) fundamental mode
### -------------------------------------------------------------------------------------------------------------------------------------------------



### CHECK PARAMETERS --------------------------------------------------------------------------------------------------------------------------------
if len(set(map(len, (soiltypes, thicknesses, Ns, fracs)))) != 1:
    raise ValueError(f"Arrays are not the same size : {soiltypes = }, {thicknesses = }, {Ns = }, {fracs = }")
  
for i, layer in enumerate(under_layers) :
    if len(layer) != 4 :
        raise ValueError(f"Under layers format is not correct:\nLayer {i+1} is '{layer}'\nCorrect format : 'thickness Vp Vs rho'")
    if i == len(under_layers) - 1 and layer[0] != 0 :
        raise ValueError(f"Last layer thickness must be 0\nLayer {i+1} is '{layer}'\nCorrect format : '0 Vp Vs rho'")
    if i < len(under_layers) - 1 and layer[0] == 0 :
        raise ValueError(f"Layers thickness must be greater than 0 (except for the last layer)\nLayer {i+1} is '{layer}'\nCorrect format : 'thickness Vp Vs rho'")
### -------------------------------------------------------------------------------------------------------------------------------------------------



# Runnung the program for each water table depth
for iWT, WT in enumerate(WTs) :
    #### ROCK PHYSICS -------------------------------------------------------------------------------------------------------------------------------
    # Saturation profile with depth
    hs, Sws, Swes = vanGen(zs, WT, soiltypes, thicknesses)

    plt.rcParams.update({'font.size': 10})
    cm = 1/2.54
    fig, axs = plt.subplots(1, 2, figsize=(8*cm, 4*cm), dpi=600, gridspec_kw={'hspace':0.5})
    axs[0].plot(hs, zs)
    axs[0].axhline(-WT, color='gray', linestyle='--')
    axs[0].set_xlabel('Pressure head h')
    axs[0].set_ylabel('Depth [m]')
    axs[0].set_ylim([zs[-1], 0])
    axs[0].grid()
    axs[0].xaxis.set_label_position('top')
    axs[0].tick_params(which='both', labelbottom=False, labeltop=True, labelleft=True, labelright=False, bottom=True, top=True, left=True, right=True)
    axs[1].plot(Sws, zs, Swes, zs)
    axs[1].axhline(-WT, color='gray', linestyle='--')
    axs[1].set_xlabel('Saturation')
    axs[1].set_ylabel('')
    axs[1].set_xlim([-0.05,1.05])
    axs[1].set_ylim([zs[-1], 0])
    axs[1].xaxis.set_label_position('top')
    axs[1].tick_params(which='both', labelbottom=False, labeltop=True, labelleft=False, labelright=False, bottom=True, top=True, left=True, right=True)
    axs[1].grid()
    fig.legend(['', 'WT', '${S_w}$', '${S_{we}}$'], loc='lower center')
    fig.savefig(f'./output/1_vanGen_WT{WT}.svg', bbox_inches='tight')

    # Effective Grain Properties (constant with depth)
    mus, ks, rhos, nus = hillsAverage(mu_clay, mu_silt, mu_sand, rho_clay,
                                      rho_silt, rho_sand, k_clay, k_silt,
                                      k_sand, soiltypes)
        
    print('\nhillsAverage')
    print(f'Shear moduli of grains: {mus = } [Pa]')
    print(f'Bulk moduli of grains: {ks = } [Pa]')
    print(f'Density of  grains: {rhos = } [kg/m3]')
    print(f"Poisson's ratio: {nus = }")
    print('\n')


    # Effective Fluid Properties
    kfs, rhofs, rhobs = effFluid(Sws, kw, ka, rhow,
                                 rhoa, rhos, soiltypes, thicknesses, dz)
    
    fig, axs = plt.subplots(1, 2, figsize=(8*cm, 4*cm), dpi=600, gridspec_kw={'hspace':0.5})
    axs[0].plot(kfs, zs)
    axs[0].axhline(-WT, color='gray', linestyle='--')
    axs[0].set_xlabel('Effective compressibility k_f [Pa-1]')
    axs[0].set_ylabel('Depth [m]')
    axs[0].set_ylim([zs[-1], 0])
    axs[0].grid()
    axs[0].xaxis.set_label_position('top')
    axs[0].tick_params(which='both', labelbottom=False, labeltop=True, labelleft=True, labelright=False, bottom=True, top=True, left=True, right=True)
    axs[1].plot(rhofs, zs, rhobs, zs)
    axs[1].axhline(-WT, color='gray', linestyle='--')
    axs[1].set_xlabel('Density rho [kg/m3]')
    axs[1].set_ylabel('')
    axs[1].set_xlim([0, 2000])
    axs[1].set_ylim([zs[-1], 0])
    axs[1].tick_params(which='both', labelbottom=False, labeltop=True, labelleft=False, labelright=False, bottom=True, top=True, left=True, right=True)
    axs[1].grid()
    axs[1].legend(['Effective fluid density rhof', 'Bulk density rhob'], fontsize=8)
    fig.savefig(f'./output/3_effFluid_WT{WT}.svg', bbox_inches='tight')


    # Hertz Mindlin Frame Properties
    KHMs, muHMs = hertzMindlin(Swes, zs, hs, rhobs,
                               g, rhoa, rhow, Ns,
                               mus, nus, fracs, kk,
                               soiltypes, thicknesses)
    
    fig, axs = plt.subplots(1, 2, figsize=(8*cm, 4*cm), dpi=600, gridspec_kw={'hspace':0.5})
    ax = axs[0]
    ax.plot(KHMs, zs)
    ax.plot(muHMs, zs)
    ax.axhline(-WT, color='gray', linestyle='--')
    ax.set_xlabel('Pressure [Pa]')
    ax.set_ylabel('Depth z [m]')
    ax.set_ylim([zs[-1], 0])
    ax.grid()
    ax.xaxis.set_label_position('top')
    ax.tick_params(which='both', labelbottom=False, labeltop=True, labelleft=True, labelright=False, bottom=True, top=True, left=True, right=True)
    ax.legend(['Effective bulk KHM', 'Shear moduli muHM'], fontsize=8)
    fig.savefig(f'./output/4_hertzMindlin_WT{WT}.svg', bbox_inches='tight')


    # Saturated Properties
    VPs, VSs = biotGassmann(KHMs, muHMs, ks, kfs,
                            rhobs, soiltypes, thicknesses, dz)
     
    fig, axs = plt.subplots(1, 2, figsize=(8*cm, 4*cm), dpi=600, gridspec_kw={'hspace':0.5})
    ax = axs[0]
    ax.plot(VPs, zs)
    ax.plot(VSs, zs)
    ax.axhline(-WT, color='gray', linestyle='--')
    ax.set_xlabel('Velocity [m/s]')
    ax.set_ylabel('Depth z [m]')
    ax.set_ylim([zs[-1], 0])
    ax.grid()
    ax.xaxis.set_label_position('top')
    ax.tick_params(which='both', labelbottom=False, labeltop=True, labelleft=True, labelright=False, bottom=True, top=True, left=True, right=True)
    ax.legend(['Vp', 'Vs'], fontsize=8)
    fig.savefig(f'./output/5_biotGassmann_WT{WT}.svg', bbox_inches='tight')
    plt.rcParams.update({'font.size': 14})


    plt.close('all')
    ### ---------------------------------------------------------------------------------------------------------------------------------------------



    #### SEISMIC FWD MODELING -----------------------------------------------------------------------------------------------------------------------
    # First arrival time computations
    thks_tmp = np.copy(thks) # Thicknesses of the layers
    VPs_tmp = np.copy(VPs) # P-wave velocities of the layers
    VSs_tmp = np.copy(VSs) # S-wave velocities of the layers
    for layer in under_layers:
        thickness = layer[0]
        if thickness == 0:
            thickness = 2*dz
        vp = layer[1]
        vs = layer[2]
        thks_tmp = np.concatenate((thks_tmp, [dz]*int(thickness/dz))) # Thicknesses of the layers
        VPs_tmp = np.concatenate((VPs_tmp, [vp]*int(thickness/dz))) # P-wave velocities of the layers
        VSs_tmp = np.concatenate((VSs_tmp, [vs]*int(thickness/dz))) # S-wave velocities of the layers            
    ThodPs = firstArrival(thks_tmp, VPs_tmp, xs, trig) # P-wave first arrival times
    ThodSs = firstArrival(thks_tmp, VSs_tmp, xs, trig) # S-wave first arrival times
    

    # Velocity model in string format for GPDC
    under_layers_str = '\n'.join([' '.join(map(str, layer)) for layer in under_layers]) + '\n'
    velocity_model_string = writeVelocityModel(thks, VPs, VSs, rhobs, under_layers_str, n_under_layers)

    # Dispersion curves computing with GPDC
    velocity_model_RAMfile = StringIO(velocity_model_string) # Keep velocity model string in the RAM in a file format alike to trick GPDC which expects a file
    gpdc_command = [f"gpdc -{wave} {n_modes} -n {nf} -min {min_f} -max {max_f} -s {s}"]

    try:
        process = run(gpdc_command, input=velocity_model_RAMfile.getvalue(), text=True, shell=True, stdout=PIPE, stderr=PIPE, check=True) # Raw output string from GPDC
    except CalledProcessError as e:
        print(f"\nERROR during GPDC computation. Returned:\n{e.stdout}")
        print("Used parameters:")
        print(f'{soiltypes = }')
        print(f'{thicknesses = }')
        print(f'{Ns = }')
        print(f'{fracs = }')
        print(f'{WT = }')
        print(f'{dz = }\n')
        print('INFO : Try to reduce dz\n')
        raise

    gpdc_output_string = process.stdout # Raw output string from GPDC
    dispersion_data, n_modes = readDispersion(gpdc_output_string) # Reads GPDC output and converts dispersion data to a list of numpy arrays for each mode
                                                                # Updates number of computed modes (can be lower than what was defined if frequency range too small)
    ### ---------------------------------------------------------------------------------------------------------------------------------------------



    ### PLOTS ---------------------------------------------------------------------------------------------------------------------------------------
    color = color_map[iWT,:] # Color for the plot


    zs_seismic = np.copy(zs)
    if n_under_layers != 0:
        for i, (thickness, Vp, Vs, rho) in enumerate(under_layers):
            thickness_number = int(max(dz, thickness) // dz)
            VSs = np.concatenate((VSs, [Vs]*thickness_number))
            VPs = np.concatenate((VPs, [Vp]*thickness_number))
            rhobs = np.concatenate((rhobs, [rho]*thickness_number))
            
            zs_to_add = np.linspace(zs_seismic[-1] - dz, zs_seismic[-1] - dz * thickness_number, thickness_number)                
            zs_seismic = np.concatenate((zs_seismic, zs_to_add))
            

    if iWT == 0:
        # Create figures
        figPage1, axs1 = plt.subplots(2, 3, figsize=(11.69,8.27), gridspec_kw={'hspace':0.3, 'wspace':0.4})
        figPage2, axs2 = plt.subplots(2, 3, figsize=(11.69,8.27), gridspec_kw={'hspace':0.3, 'wspace':0.4})

        # Name for the output files with soil types and thicknesses
        name = ''
        for soiltype, thickness in zip(soiltypes, thicknesses):
            name += f"{soiltype}{thickness}_"
        name = name[:-1]
        if n_under_layers != 0:
            name += "_substratum"



    # Plot VPs vs. zs
    axs1[0,0].plot(VPs, zs_seismic, linewidth=1.5, color=color)
    axs1[0,0].axhline(-WT, color=color, linestyle='--', linewidth=0.5)
    axs1[0,0].set_xlim(0, max(VPs) + 100)
    axs1[0,0].set_xlabel('$V_p$ [m/s]')
    axs1[0,0].set_ylabel('$z$ [m]')

    # Plot VSs vs. zs
    axs1[0,1].plot(VSs, zs_seismic, linewidth=1.5, color=color)
    axs1[0,1].axhline(-WT, color=color, linestyle='--', linewidth=0.5)
    axs1[0,1].set_xlim(0, max(VSs) + 100)
    axs1[0,1].set_xlabel('$V_s$ [m/s]')
    axs1[0,1].set_ylabel('$z$ [m]')

    # Plot rho_b vs. zs
    axs1[0,2].plot(rhobs, zs_seismic, linewidth=1.5, color=color)
    axs1[0,2].axhline(-WT, color=color, linestyle='--', linewidth=0.5)
    axs1[0,2].set_xlim(min(rhobs) - 100, max(rhobs) + 100)
    axs1[0,2].set_xlabel('$\\rho_b$ [kg/$m^3$]')
    axs1[0,2].set_ylabel('$z$ [m]')

    # Plot Sws vs. zs
    axs1[1,0].plot(Sws, zs, linewidth=1.5, color=color)
    axs1[1,0].axhline(-WT, color=color, linestyle='--', linewidth=0.5)
    axs1[1,0].set_xlim(0, 1.1)
    axs1[1,0].set_xlabel('$S_w$')
    axs1[1,0].set_ylabel('$z$ [m]')

    # Plot Poisson ratio vs. zs
    poisson_ratios = list(map(fish, VPs, VSs))
    axs1[1,1].plot(poisson_ratios, zs_seismic, linewidth=1.5, color=color)
    axs1[1,1].axhline(-WT, color=color, linestyle='--', linewidth=0.5)
    axs1[1,1].set_xlim(0, max(poisson_ratios) + 0.1)
    axs1[1,1].set_xlabel('Poisson ratio')
    axs1[1,1].set_ylabel('$z$ [m]')

    # No plot
    axs1[1,2].axis('off')



    # Plot VPs vs. zs
    axs2[0,0].plot(VPs, zs_seismic, linewidth=1.5, color=color)
    axs2[0,0].axhline(-WT, color=color, linestyle='--', linewidth=0.5)
    axs2[0,0].set_xlim(0, max(VPs) + 100)
    axs2[0,0].set_xlabel('$V_p$ [m/s]')
    axs2[0,0].set_ylabel('$z$ [m]')

    # Plot VSs vs. zs
    axs2[0,1].plot(VSs, zs_seismic, linewidth=1.5, color=color)
    axs2[0,1].axhline(-WT, color=color, linestyle='--', linewidth=0.5)
    axs2[0,1].set_xlim(0, max(VSs) + 100)
    axs2[0,1].set_xlabel('$V_s$ [m/s]')
    axs2[0,1].set_ylabel('$z$ [m]')

    # Plot rho_b vs. zs
    axs2[0,2].plot(rhobs, zs_seismic, linewidth=1.5, color=color)
    axs2[0,2].axhline(-WT, color=color, linestyle='--', linewidth=0.5)
    axs2[0,2].set_xlim(min(rhobs) - 100, max(rhobs) + 100)
    axs2[0,2].set_xlabel('$\\rho_b$ [kg/$m^3$]')
    axs2[0,2].set_ylabel('$z$ [m]')

    # Plot simulated P-wave first arrivals
    axs2[1,0].plot(xs, ThodPs, linewidth=1.5, color=color)
    axs2[1,0].set_xlim([0, max(xs)])
    axs2[1,0].set_ylim([0, max(ThodPs)])
    axs2[1,0].set_xlabel('Offset [m]')
    axs2[1,0].set_ylabel('P- first arrival time [s]')

    # Plot simulated S-wave first arrivals
    axs2[1,1].plot(xs, ThodSs, linewidth=1.5, color=color)
    axs2[1,1].set_xlim([0, max(xs)])
    axs2[1,1].set_ylim([0, max(ThodSs)])
    axs2[1,1].set_xlabel('Offset [m]')
    axs2[1,1].set_ylabel('S- first arrival time [s]')

    # Plot simulated dispersion curve
    max_vr = np.max(dispersion_data[0][:,1])
    min_vr = np.min(dispersion_data[0][:,1])
    for mode in range(n_modes):
        axs2[1,2].plot(dispersion_data[mode][:,0], dispersion_data[mode][:,1], linewidth=1.5, color=color)
        if np.max(dispersion_data[mode][:,1]) > max_vr:
            max_vr = np.max(dispersion_data[mode][:,1])
        if np.min(dispersion_data[mode][:,1]) < min_vr:
            min_vr = np.min(dispersion_data[mode][:,1])
        # np.savetxt(f'./output/{name}_M{mode}_WT{WT}.txt', dispersion_data[mode], fmt='%f') # Save dispersion data to file
    axs2[1,2].set_xlim([min_f-5, max_f+5])
    # axs2[1,2].set_ylim([min_vr-100, max_vr+100])
    axs2[1,2].set_ylim([100, 500])
    axs2[1,2].set_xlabel('Frequency [Hz]')
    axs2[1,2].set_ylabel('P-SV phase vel. [m/s]')
    ### ---------------------------------------------------------------------------------------------------------------------------------------------



### CONVERT PLOTS TO SINGLE PDF ---------------------------------------------------------------------------------------------------------------------
pdf_filename = f"./output/{name}.SANTILUDO.pdf"
pdf_pages = PdfPages(pdf_filename)

pdf_pages.savefig(figPage1, bbox_inches='tight')
pdf_pages.savefig(figPage2, bbox_inches='tight')

pdf_pages.close()
### -------------------------------------------------------------------------------------------------------------------------------------------------



### TIME END ----------------------------------------------------------------------------------------------------------------------------------------
end = perf_counter()
print(f"\nElapsed time : {end-start} s\n")
### -------------------------------------------------------------------------------------------------------------------------------------------------

