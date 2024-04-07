# PROGRAM SANTILUDO*
- Programmers: *S.G. Solazzi & L. Bodet* (function *Prof_VanGenLB.m* courtesy of D. Jougnot)
- Frist version: *2020/07/14*
- Last update (first draft with synthetic computation): *2022/03*
- some minor modifs 2024/03
- Converted to C++, and interfaced with Python, on 2024/04 by J. Cunha Teixeira

## Objectives:
This program takes the hydromechanical properties of a partially saturated 1D soil for a given steady-state condition and computes:
- (i)   The saturation profile in depth ;
- (ii)  The variations of Vp, Vs, and rho with saturation (depth) ;    
- (iii) The P-wave first arrival times and suface-wave dispersion for a given linear acquisition setup.

## Reference
Solazzi, S. G., Bodet, L., Holliger, K., & Jougnot, D. (2021). Surface-wave dispersion in partially saturated soils: The role of capillary forces. Journal of Geophysical Research: Solid Earth, 126, e2021JB022074. 

https://doi.org/10.1029/2021JB022074

-> this code makes it possible to compute figures 5 and 6 (but also the others, if you manage some minor modifs...)

## Dependencies:
  Source code is stored in ./src folder.
  Functions are written in C++ file (*.cpp) and wrapped in Cython file (*.pyx) to be compiled
  as a shared library and called in Python.
  - Source functions in C++ are named as nameFunction_src
  - Wapped functions in Cython to be called in Python are named as nameFunction
  -> If you change the header of a function in the *.cpp, you must change it on the *.pyx file also.

  ### C++ functions for RP and VG models are stored in ./src/RPfunctions_src.cpp and ./src/GVfunctions_src.cpp:
      - ./src/GVfunctions_src.cpp/vanGen_src -> to compute Saturation from...
      - ./src/GVfunctions_src.cpp/hillsAverage_src -> to compute ... from ...
      - ...
  ### Seismic data forward modelling utils stored in ./src/TTDSPfunctions_src.cpp :
      - ./src/TTDSPfunctions_src.cpp/firstArrival_src -> compute P- or S-wave first arrival times from a 1D velocity model
      + *gpdc* -> code (to be compiled and called) from https://www.geopsy.org 
          to to compute surface-wave dispersion from a 1D velocity model
      + I/O functions:
      - ./src/TTDSPfunctions_src.cpp/dinSave_src -> create 1D velocity models readable by *gpdc*
      - ./src/TTDSPfunctions_src.cpp/readDisp -> Only function not wtitten in C++ but directly in Cython in the *.pyx file - reads dispersion curves created by *gpdc*
 ### Utils in bash and or latex to plot results ./src folder:
      - *cropSANTILUDO.sh* (as to be chmod +x...) 
      - ...

## Compilation is done with the setup.py file :
python3 setup.py build_ext --inplace clean
  
This creates a ./build and ./bib directories, with the last containing the shared libraries *.so
that can be used in Python and the resulting *.cpp files from the *.pyx compilation

## Running the main script :
python3 main_SANTILUDO.py

## Comments:
This is a quickly written prototype. A lot to do to improve this code :
- add/improve comments (names, variables, units, inputs/outputs, 
  required function and called codes, dependencies etc)
- add references to associated litterature
- optimise/vectorise if needed
- transform repetitive actions into functions !!!

## Identified problems:
- less Cells are needed when soiltype is "sand"
- the definition of variable *xs* versus variable *Xobs* is unclear
- the change of depth scale is weird
- the mixture proportions and associated RP param is to be clarified 
- ...

## 2DO list: 
- we have to test the sensitivity to layer thikness, half-space depth and 
  thickness, layer discretisation *etc* !!!!!
- we have to test the sensitivity to HM parameters !!!
- we have to test the sensitivity to meca parameters and composition !!!
- we have to test the sensitivity geom parameters !!!

- we have to implement and alternative RP models !!!!!
- ...

## About us
<kbd>Under construction...</kbd>

-- 
Dr Santiago G. Solazzi<br />
Researcher (Subsoil Technologies)<br />
YPF-Tecnología, Argentina<br />
Verified email at ypftecnologia.com<br />
https://orcid.org/0000-0002-6952-0309<br />

-- 
Dr Ludovic Bodet (ludovic.bodet sorbonne-universite.fr) <br />
Associate Professor/Maître de conférences, HDR [https://hal.sorbonne-universite.fr/tel-02866882v2](https://hal.sorbonne-universite.fr/tel-02866882v2]) <br />
Head of the Hydrology, Hydrogeology and Near-Surface Geophysics' department at sorbonne-universite.fr <br />
Head of the Critical Zone Hydro&Geophysical Instrument Toolbox (Terre-Mer-Sols/OSU ECCETERRA/SU) <br />

UMR CNRS METIS <br />
Sorbonne Université, Tours 46-56, étage 3, case courrier 105 <br />
4 place Jussieu, 75252 Paris cedex 05, FRANCE <br />
