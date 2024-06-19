# cython: language_level=3
# distutils: language=c++
# cython: c_string_type=unicode, c_string_encoding=utf8

# Import necessary Cython declarations
from libcpp.vector cimport vector
from libcpp.string cimport string


# Import the C++ function declaration
cdef extern from "RPfunctions_src.cpp":

    cdef struct effFluidResult:
       vector[double] k_f, rho_f, rho_b

    cdef struct biotGassmannResult:
        vector[double] v_p, v_s

    cdef struct hertzMindlinResult:
        vector[double] k_m, mu_m

    cdef struct hillsAverageResult:
        vector[double] mu_s, k_s, rho_s, nu_s

    effFluidResult effFluidSrc(vector[double] s_w, double k_w, double k_a, double rho_w, double rho_a, vector[double] rho_s, vector[string] soil_types, vector[double] thicknesses, double dz)
    double fishSrc(double v_p, double v_s)
    biotGassmannResult biotGassmannSrc(vector[double] k_m, vector[double] mu_m, vector[double] k_s, vector[double] k_f, vector[double] rho_b, vector[string] soil_types, vector[double] thicknesses, double dz)
    hertzMindlinResult hertzMindlinSrc(vector[double] s_we, vector[double] z, vector[double] h, vector[double] rho_b, double g, double rho_a, double rho_w, vector[double] n_s, vector[double] mu_s, vector[double] nu_s, vector[double] fracs, int kk, vector[string] soil_types, vector[double] thicknesses)
    hillsAverageResult hillsAverageSrc(double mu_clay, double mu_silt, double mu_sand, double rho_clay, double rho_silt, double rho_sand, double k_clay, double k_silt, double k_sand, vector[string] soil_types)



# Define Python wrapper functions
def effFluid(s_w, k_w, k_a, rho_w, rho_a, rho_s, soil_types, thicknesses, dz):
    cdef effFluidResult result = effFluidSrc(s_w, k_w, k_a, rho_w, rho_a, rho_s, soil_types, thicknesses, dz)
    return (result.k_f, result.rho_f, result.rho_b)

def fish(v_p, v_s):
    cdef double sig = fishSrc(v_p, v_s)
    return sig

def biotGassmann(k_m, mu_m, k_s, k_f, rho_b, soil_types, thicknesses, dz):
    cdef biotGassmannResult result = biotGassmannSrc(k_m, mu_m, k_s, k_f, rho_b, soil_types, thicknesses, dz)
    return (result.v_p, result.v_s)

def hertzMindlin(s_we, z, h, rho_b, g, rho_a, rho_w, n_s, mu_s, nu_s, fracs, kk, soil_types, thicknesses):
    cdef hertzMindlinResult result = hertzMindlinSrc(s_we, z, h, rho_b, g, rho_a, rho_w, n_s, mu_s, nu_s, fracs, kk, soil_types, thicknesses)
    return (result.k_m, result.mu_m)

def hillsAverage(mu_clay, mu_silt, mu_sand, rho_clay, rho_silt, rho_sand, k_clay, k_silt, k_sand, soil_types):
    cdef hillsAverageResult result = hillsAverageSrc(mu_clay, mu_silt, mu_sand, rho_clay, rho_silt, rho_sand, k_clay, k_silt, k_sand, soil_types)
    return (result.mu_s, result.k_s, result.rho_s, result.nu_s)