# cython: language_level=3
# distutils: language=c++
# cython: c_string_type=unicode, c_string_encoding=utf8

# Import necessary Cython declarations
from libcpp.vector cimport vector
from libcpp.string cimport string



# Import the C++ function declaration
cdef extern from "VGfunctions_src.cpp":

    struct vanGenResult:
        vector[double] h, s_w, s_we

    struct soilType:
        double wsand, wclay, wsilt, phi, alpha, nvg, theta, s_wr

    vanGenResult vanGenSrc(vector[double] z, double wt, vector[string] soil_types, vector[double] thicknesses)
    soilType selectSoilTypeSrc(string soil_type)



# Define Python wrapper functions
def vanGen(vector[double] z, double wt, vector[string] soil_types, vector[double] thicknesses):
    cdef vanGenResult result = vanGenSrc(z, wt, soil_types, thicknesses)
    return result.h, result.s_w, result.s_we

def selectSoilType(string soil_type):
    cdef soilType result = selectSoilTypeSrc(soil_type)
    return result.wsand, result.wclay, result.wsilt, result.phi, result.alpha, result.nvg, result.theta, result.s_wr 
