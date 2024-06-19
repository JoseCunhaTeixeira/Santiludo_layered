#include <cmath>
#include <vector>
#include "VGfunctions_src.h"



struct effFluidResult {
    std::vector<double> k_f; // Effective compressibility over depth [Pa^-1]
    std::vector<double> rho_f; // Effective fluid density over depth [Kg/m3]
    std::vector<double> rho_b; // Bulk density over depth [Kg/m3]
};

struct biotGassmannResult {
    std::vector<double> v_p; // P-wave velocity [m/s]
    std::vector<double> v_s; // S-wave velocity [m/s]
};

struct hertzMindlinResult {
    std::vector<double> k_m; // Effective bulk modulus over depth [Pa]
    std::vector<double> mu_m; // Effective shear modulus over depth [Pa]
};

struct hillsAverageResult {
    std::vector<double> mu_s; // Shear moduli of grain of each layer [Pa]
    std::vector<double> k_s; // Bulk moduli of grain of each layer [Pa]
    std::vector<double> rho_s; // Densities of grain of each layer [Kg/m3]
    std::vector<double> nu_s; // Poisson's ratios of each layer [-]
};


/**
 * Computes the effective properties of the solid grains from its constituents for each layer of soil types
 * 
 * @author S.G. Solazzi
 * @author J. Cunha Teixeira
 * @param mu_clay shear modulus of clay [GPa]
 * @param mu_silt shear modulus of silt [GPa]
 * @param mu_sand shear modulus of sand [GPa]
 * @param rho_clay density of clay [Kg/m3]
 * @param rho_silt density of silt [Kg/m3]
 * @param rho_sand density of sand [Kg/m3]
 * @param k_clay bulk modulus of clay [GPa]
 * @param k_silt bulk modulus of silt [GPa]
 * @param k_sand bulk modulus of sand [GPa]
 * @param soil_types soil types of each layer in the profile
 * @return a hillsAverageResult structure containing the effective properties of the solid grains for each layer of soil types
*/
// 2021 - Solazzi - Implemented in Matlab
// MAR2024 - Cunha Teixeira - Traduced in C++
// APR2024 - Cunha Teixeira - Adapted for mutli-layers
hillsAverageResult hillsAverageSrc(const double& mu_clay,
                                     const double& mu_silt,
                                     const double& mu_sand,
                                     const double& rho_clay,
                                     const double& rho_silt,
                                     const double& rho_sand,
                                     const double& k_clay,
                                     const double& k_silt,
                                     const double& k_sand,
                                     const std::vector<std::string>& soil_types) {

    std::vector<double> mu_s(soil_types.size()); // Shear moduli of grain of each layer
    std::vector<double> k_s(soil_types.size()); // Bulk moduli of grain of each layer
    std::vector<double> rho_s(soil_types.size()); // Densities of grain of each layer
    std::vector<double> nu_s(soil_types.size()); // Poisson's ratios of each layer

    for (size_t j = 0; j < soil_types.size(); ++j) {
        std::string_view soil_type = soil_types[j]; // Soil type for current layer

        soilType soil = selectSoilTypeSrc(soil_type); // Soil properties for current layer

        // Shear moduli of grain [Pa]
        mu_s[j] = 0.5 * (1 / (soil.wclay / mu_clay + soil.wsilt / mu_silt + soil.wsand / mu_sand) + (soil.wclay * mu_clay + soil.wsilt * mu_silt + soil.wsand * mu_sand)) * 1e9;
        
        // Bulk moduli of grain [Pa]
        k_s[j] = 0.5 * (1 / (soil.wclay / k_clay + soil.wsilt / k_silt + soil.wsand / k_sand) + (soil.wclay * k_clay + soil.wsilt * k_silt + soil.wsand * k_sand)) * 1e9;
        
        // Densities of grain [Kg/m3]
        rho_s[j] = soil.wclay * rho_clay + soil.wsilt * rho_silt + soil.wsand * rho_sand;
        
        // Poisson's ratios
        nu_s[j] = (3 * k_s[j] - 2 * mu_s[j]) / (2 * (3 * k_s[j] + mu_s[j]));
    }

    hillsAverageResult result; // Structure to store the results
    result.mu_s = mu_s;
    result.k_s = k_s;
    result.rho_s = rho_s;
    result.nu_s = nu_s;

    return result;
}


/**
 * Computes effective fluid properties and bulk density for each layer of soil type
 * 
 * @author S.G. Solazzi
 * @author J. Cunha Teixeira
 * @param s_ws total saturation over depth [-]
 * @param kw water bulk modulus [Pa]
 * @param ka air bulk modulus [Pa]
 * @param rho_w water density [Kg/m3]
 * @param rho_a air density [Kg/m3]
 * @param rho_s Densities of grain of each layer [Kg/m3]
 * @param soil_types Soil types of each layer in the profile
 * @param thiknesses thickness of each layer [m]
 * @param dz depth discretization [m]
 * @return a effFluidResult structure containing the effective properties and bulk density for each layer of soil type
*/
// 2021 - Solazzi - Implemented in Matlab
// MAR2024 - Cunha Teixeira - Traduced in C++
// APR2024 - Cunha Teixeira - Adapted for mutli-layers
effFluidResult effFluidSrc(const std::vector<double>& s_ws,
                             const double& kw,
                             const double& ka,
                             const double& rho_w,
                             const double& rho_a,
                             const std::vector<double>& rho_s,
                             const std::vector<std::string>& soil_types,
                             const std::vector<double>& thicknesses,
                             const double& dz) {

    size_t num_points = s_ws.size(); // Number of cells in the profile

    std::vector<double> k_f(num_points); // Effective compressibility over depth
    std::vector<double> rho_f(num_points); // Effective fluid density over depth
    std::vector<double> rho_b(num_points); // Bulk density over depth

    int start = 0; // Depth start index for fisrt layer

    for (size_t j = 0; j < soil_types.size(); ++j) {
        std::string_view soil_type = soil_types[j]; // Soil type for current layer
        double thickness = thicknesses[j]; // Thickness of current layer

        soilType soil = selectSoilTypeSrc(soil_type); // Soil properties for current layer

        int end = round(thickness / dz) + start; // Depth end index for current layer

        for (int i = start; i < end; ++i) {
            k_f[i] = 1.0 / (s_ws[i] / kw + (1.0 - s_ws[i]) / ka); // Effective compressibility (Woods formula)
            rho_f[i] = s_ws[i] * rho_w + (1.0 - s_ws[i]) * rho_a; // Effective fluid density
            rho_b[i] = (1.0 - soil.phi) * rho_s[j] + soil.phi * rho_f[i]; // Bulk density
        }
        start = end; // Update start index for next layer
    }

    effFluidResult result; // Structure to store the results
    result.k_f = k_f;
    result.rho_f = rho_f;
    result.rho_b = rho_b;

    return result;
}


/**
 * Computes the effective properties for each layer of soil type using the Hertz Mindlin model
 * 
 * @author S.G. Solazzi
 * @author J. Cunha Teixeira
 * @param s_we effective wetting phase saturation over depth [-]
 * @param z depths of the profile [m]
 * @param h pressure head over depth [m]
 * @param rho_b bulk density over depth [Kg/m3]
 * @param g gravity [m/s2]
 * @param rho_a air density [Kg/m3]
 * @param rho_w water density [Kg/m3]
 * @param n_s coordination number of each layer [-]
 * @param mu_s shear moduli of grain of each layer [GPa]
 * @param nu_s Poisson's ratios of each layer [-]
 * @param fracs fraction of non-slipping grains of each layer [-]
 * @param kk type of effective stress
 * @param soil_types soil types of each layer in the profile
 * @param thicknesses thickness of each layer [m]
 * @return a hertzMindlinResult structure containing the effective properties of each layer of soil type
*/
// 2021 - Solazzi - Implemented in Matlab
// MAR2024 - Cunha Teixeira - Traduced in C++
// APR2024 - Cunha Teixeira - Adapted for mutli-layers
hertzMindlinResult hertzMindlinSrc(const std::vector<double>& s_we,
                                     const std::vector<double>& z,
                                     const std::vector<double>& h,
                                     const std::vector<double>& rho_b,
                                     const double& g,
                                     const double& rho_a,
                                     const double& rho_w,
                                     const std::vector<double>& n_s,
                                     const std::vector<double>& mu_s,
                                     const std::vector<double>& nu_s,
                                     const std::vector<double>& fracs,
                                     const int& kk,
                                     const std::vector<std::string>& soil_types,
                                     const std::vector<double>& thicknesses) {

    size_t num_points = s_we.size(); // Number of cells in the profile

    std::vector<double> k_m(num_points); // Effective bulk modulus over depth
    std::vector<double> mu_m(num_points); // Effective shear modulus over depth
    
    std::vector<double> pe(num_points); // Overburden stress over depth

    double dz = fabs(z[1] - z[0]); // Depth discretization

    int start = 0; // Depth start index for fisrt layer

    double sigma_prev = 0; // Previous overburden stress

    for (size_t j = 0; j < soil_types.size(); ++j) {
        std::string_view soil_type = soil_types[j]; // Soil type for current layer
        double thickness = thicknesses[j]; // Thickness of current layer

        soilType soil = selectSoilTypeSrc(soil_type); // Soil properties for current layer

        int end = round(thickness / dz) + start; // Depth end index for current layer

        for (int i = start; i < end; ++i) {

            double sigma = rho_b[i] * g * dz + sigma_prev; // Overburden stress
            sigma_prev = sigma;

            double dep = fabs(z[i]); // Depth (positive value)
            double st = (h[i] <= 0) ? 0 : s_we[i] * rho_w * g * h[i]; // Suction term | Null when fully saturated
            double pa = rho_a * g * dep; // Air pressure
            double pe; // Effective stress

            if (kk == 1) {
                pe = 101325;
            } else if (kk == 2) {
                pe = sigma - pa;
            } else if (kk == 3) {
                pe = sigma - pa + st;
            } else {
                throw std::invalid_argument("Invalid value for kk");
            }

            k_m[i] = std::pow((n_s[j]*n_s[j] * std::pow((1 - soil.phi), 2) * mu_s[j]*mu_s[j] * pe) / (18 * M_PI*M_PI * std::pow((1 - nu_s[j]), 2)), 1.0 / 3.0); // Effective bulk modulus
            mu_m[i] = ((2 + 3 * fracs[j] - (1 + 3 * fracs[j]) * nu_s[j]) / (5 * (2 - nu_s[j]))) * std::pow((3 * n_s[j]*n_s[j] * std::pow((1 - soil.phi), 2) * mu_s[j]*mu_s[j] * pe) / (2 * M_PI*M_PI * std::pow((1 - nu_s[j]), 2)), 1.0 / 3.0); // Effective shear modulus
        }
        start = end; // Update start index for next layer
    }

    hertzMindlinResult result; // Structure to store the results
    result.k_m = k_m;
    result.mu_m = mu_m;
    
    return result;
}


/**
 * Computes the P- and S-wave velocity for each layer of soil type using the Biot Gassmann model
 * 
 * @author S.G. Solazzi
 * @author J. Cunha Teixeira
 * @param k_m effective bulk modulus over depth [Pa]
 * @param mu_m effective shear modulus over depth [Pa]
 * @param k_s bulk moduli of grain of each layer [Pa]
 * @param k_f effective compressibility over depth [Pa^-1]
 * @param rho_b bulk density over depth [Kg/m3]
 * @param soil_types soil types of each layer in the profile
 * @param thicknesses thickness of each layer [m]
 * @param dz depth discretization [m]
 * @return a biotGassmannResult structure containing the effective properties of each layer of soil type
*/
// 2021 - Solazzi - Implemented in Matlab
// MAR2024 - Cunha Teixeira - Traduced in C++
// APR2024 - Cunha Teixeira - Adapted for mutli-layers
biotGassmannResult biotGassmannSrc(const std::vector<double>& k_m,
                                     const std::vector<double>& mu_m,
                                     const std::vector<double>& k_s,
                                     const std::vector<double>& k_f,
                                     const std::vector<double>& rho_b,
                                     const std::vector<std::string>& soil_types,
                                     const std::vector<double>& thicknesses,
                                     const double& dz) {

    size_t num_points = k_m.size(); // Number of cells in the profile

    std::vector<double> v_p(num_points); // P-wave velocity over depth
    std::vector<double> v_s(num_points); // S-wave velocity over depth

    int start = 0; // Depth start index for fisrt layer

    for (size_t j = 0; j < soil_types.size(); ++j) {
        std::string_view soil_type = soil_types[j]; // Soil type for current layer
        double thickness = thicknesses[j]; // Thickness of current layer

        soilType soil = selectSoilTypeSrc(soil_type); // Soil properties for current layer

        int end = round(thickness / dz) + start; // Depth end index for current layer

        for (int i = start; i < end; ++i) {
            double alpha2, k;
            alpha2 = 1 - k_m[i] / k_s[j];
            k = k_m[i] + (alpha2 * alpha2) / (soil.phi / k_f[i] + (1 - soil.phi) / k_s[j] - k_m[i] / (k_s[j] * k_s[j]));

            v_p[i] = sqrt((k + (4.0 / 3.0) * mu_m[i]) / rho_b[i]); // P-wave velocity
            v_s[i] = sqrt(mu_m[i] / rho_b[i]); // S-wave velocity
        }
        start = end; // Update start index for next layer
    }

    biotGassmannResult result; // Structure to store the results
    result.v_p = v_p;
    result.v_s = v_s;

    return result;
}


/**
 * Computes Poisson's ratio from v_p and v_s
 * 
 * @author S.G. Solazzi
 * @author J. Cunha Teixeira
 * @param v_p P-wave velocity [m/s]
 * @param v_s S-wave velocity [m/s]
 * @return the Poisson's ratio from v_p and v_s
*/
// 2021 - Solazzi - Implemented in Matlab
// MAR2024 - Cunha Teixeira - Traduced in C++
// APR2024 - Cunha Teixeira - Adapted for mutli-layers
double fishSrc(double v_p, double v_s) {
    double ratio = v_p/v_s;
    // return (v_p*v_p - 2*v_s*v_s) / (2 * (v_p*v_p - v_s*v_s))
    return (0.5 * ratio*ratio - 1) / (ratio*ratio - 1);
}