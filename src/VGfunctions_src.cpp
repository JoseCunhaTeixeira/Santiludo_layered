#include <cmath>
#include <stdexcept>
#include <vector>
#include <string>



struct vanGenResult{
    std::vector<double> h; // pressure head over depth [m]
    std::vector<double> s_w; // total saturation over depth [-]
    std::vector<double> s_we; // effective wetting phase saturation over depth [-]
};

struct soilType{
    double wsand; // sand ratio [-]
    double wclay; // clay ratio [-]
    double wsilt; // silt ratio [-]
    double phi; // soil porosity [-]
    double alpha; // inverse of entry pressure [1/m]
    double nvg; // Van Genuchten parameter [-]
    double theta; // Van Genuchten parameter [-]
    double s_wr; // residual water saturation [-]
};



/**
 * Set the correct parameters values to a soilType for a given soil_type name
 * Based on:
 * - Carsel, R. F., & Parrish, R. S. (1988). 
 *      Developing joint probability distributions of soil water retention 
 *      characteristics. Water Resource Research, 24(5), 755–769.
 *      https://doi.org/10.1029/wr024i005p00755
 * 
 * @author S.G. Solazzi
 * @author J. Cunha Teixeir
 * @param soil_type soil type name
 * @return a selectSoilTypeResult structure containing the properties of the given soil type
*/
// 2021 - Solazzi - Implemented in Matlab
// MAR2024 - Cunha Teixeira - Traduced in C++
soilType selectSoilTypeSrc(const std::string_view& soil_type) {
    soilType soil; // Structure to store the results

    if (soil_type == "clay") {
        soil.wsand = 0.149;
        soil.wclay = 0.552;
        soil.wsilt = 1 - soil.wsand - soil.wclay;
        soil.phi = 0.38;
        soil.alpha = 0.8;
        soil.nvg = 1.09;
        soil.theta = 0.068;
        soil.s_wr = soil.theta / soil.phi;    }
    else if (soil_type == "silt") {
        soil.wsand = 0.058;
        soil.wclay = 0.095;
        soil.wsilt = 1 - soil.wsand - soil.wclay;
        soil.phi = 0.46;
        soil.alpha = 1.6;
        soil.nvg = 1.37;
        soil.theta = 0.034;
        soil.s_wr = soil.theta / soil.phi;    }
    else if (soil_type == "clayloam") {
        soil.wsand = 0.298;
        soil.wclay = 0.326;
        soil.wsilt = 1 - soil.wsand - soil.wclay;
        soil.phi = 0.41;
        soil.alpha = 1.9;
        soil.nvg = 1.31;
        soil.theta = 0.095;
        soil.s_wr = soil.theta / soil.phi;    }
    else if (soil_type == "loam") {
        soil.wsand = 0.4;
        soil.wclay = 0.197;
        soil.wsilt = 1 - soil.wsand - soil.wclay;
        soil.phi = 0.43;
        soil.alpha = 3.6;
        soil.nvg = 1.56;
        soil.theta = 0.078;
        soil.s_wr = soil.theta / soil.phi;    }
    else if (soil_type == "loamysand") {
        soil.wsand = 0.809;
        soil.wclay = 0.064;
        soil.wsilt = 1 - soil.wsand - soil.wclay;
        soil.phi = 0.41;
        soil.alpha = 12.4;
        soil.nvg = 1.28;
        soil.theta = 0.057;
        soil.s_wr = soil.theta / soil.phi;    }
    else if (soil_type == "cleansand") {
        soil.wsand = 1;
        soil.wclay = 0;
        soil.wsilt = 0;
        soil.phi = 0.43;
        soil.alpha = 14.5;
        soil.nvg = 2.68;
        soil.theta = 0.045;
        soil.s_wr = soil.theta / soil.phi;    }
    else if (soil_type == "sand") {
        soil.wsand = 0.927;
        soil.wclay = 0.029;
        soil.wsilt = 1 - soil.wsand - soil.wclay;
        soil.phi = 0.43;
        soil.alpha = 14.5;
        soil.nvg = 2.68;
        soil.theta = 0.045;
        soil.s_wr = soil.theta / soil.phi;    }
    else if (soil_type == "sandyclay") {
        soil.wsand = 0.52;
        soil.wclay = 0.43;
        soil.wsilt = 0.05;
        soil.phi = 0.38;
        soil.alpha = 2.7;
        soil.nvg = 1.23;
        soil.theta = 0.1;
        soil.s_wr = soil.theta / soil.phi;    }
    else if (soil_type == "sandyclayloam") {
        soil.wsand = 0.543;
        soil.wclay = 0.274;
        soil.wsilt = 1 - soil.wsand - soil.wclay;
        soil.phi = 0.43;
        soil.alpha = 1.0;
        soil.nvg = 1.23;
        soil.theta = 0.089;
        soil.s_wr = soil.theta / soil.phi;    }
    else if (soil_type == "sandyloam") {
        soil.wsand = 0.634;
        soil.wclay = 0.111;
        soil.wsilt = 1 - soil.wsand - soil.wclay;
        soil.phi = 0.41;
        soil.alpha = 7.5;
        soil.nvg = 1.89;
        soil.theta = 0.065;
        soil.s_wr = soil.theta / soil.phi;    }
    else {
        throw std::invalid_argument("Invalid soiltype");
    }

    return soil;
}


/**
 * Computes for each layer of soil type the properties of the soil using the van Genuchten model
 * 
 * @author S.G. Solazzi
 * @author J. Cunha Teixeir
 * @param z depths of the profile [m]
 * @param wt depth of the water table [m]
 * @param soil_types soil types of each layer in the profile
 * @param thicknesses thickness of each layer [m]
 * @return a vanGenResult structure containing the properties of each layer of soil type
*/
// 2021 - Solazzi - Implemented in Matlab
// MAR2024 - Cunha Teixeira - Traduced in C++
// APR2024 - Cunha Teixeira - Adapted for mutli-layers
vanGenResult vanGenSrc(const std::vector<double>& z,
                         const double& wt,
                         const std::vector<std::string>& soil_types,
                         const std::vector<double>& thicknesses) {

    size_t numPoints = z.size(); // Number of cells in the profile
    double dz = std::abs(z[1] - z[0]); // Depth discretization

    std::vector<double> h(numPoints); // Pressure head over depth
    std::vector<double> s_we(numPoints); // Effective wetting phase saturation over depth
    std::vector<double> s_w(numPoints); // Total saturation over depth

    int start = 0; // Depth start index for fisrt layer

    for (size_t j = 0; j < soil_types.size(); ++j) {
        std::string_view soil_type = soil_types[j]; // Soil type for current layer
        double thickness = thicknesses[j]; // Thickness of current layer
        
        int end = round(thickness / dz) + start; // Depth end index for current layer

        soilType soil = selectSoilTypeSrc(soil_type); // Soil properties for current layer
        double m = 1.0 - 1.0/soil.nvg; // Parameters related to the pore size distribution (Carsel & Parrish, 1988)

        for (int i = start; i < end; ++i) {
            h[i] = z[i] + wt;
            s_we[i] = (h[i] <= 0) ? 1 : std::pow(1 + std::pow(soil.alpha * h[i], soil.nvg), -m);
            s_w[i] = s_we[i] * (1 - soil.s_wr) + soil.s_wr;
        }
        start = end; // Update start index for next layer
    }

    vanGenResult result; // Structure to store the results
    result.h = h;
    result.s_w = s_w;
    result.s_we = s_we;

    return result;
}