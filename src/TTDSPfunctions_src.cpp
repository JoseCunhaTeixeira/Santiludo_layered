#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>


/**
 * Returns a string of the velocity model in a format expected by GPDC
 * 
 * @author S.G. Solazzi
 * @author J. Cunha Teixeira
 * @param thk thickness of each layer [m]
 * @param v_p P-wave velocity of each layer [m/s]
 * @param v_s S-wave velocity of each layer [m/s]
 * @param rho density of each layer [kg/m^3]
 * @param under_layers layers to put under the studied soil column on the velocity model in GPDC format
 * @param n_under_layers number of layers in the under_layers
 * @return a string of the velocity model in a format expected by GPDC
*/
// 2021 - Solazzi - Implemented in Matlab
// MAR2024 - Cunha Teixeira - Traduced in C++
std::string writeVelocityModelSrc(const std::vector<double>& thk,
                        const std::vector<double>& v_p,
                        const std::vector<double>& v_s,
                        const std::vector<double>& rho,
                        const std::string_view& under_layers,
                        const int& n_under_layers) {

    int nl = thk.size(); // Number of layers
    std::stringstream velocity_model; // String to store the velocity model
    velocity_model << nl + n_under_layers << '\n'; // Number of layers in the velocity model

    // Write thickness and velocity in dinver format
    for (int i = 0; i < nl; ++i) {
        if (i == nl - 1 && under_layers.size() == 0) {
            velocity_model << "0 ";
        } else {
            velocity_model << thk[i] << ' ';
        }
        velocity_model << v_p[i] << ' ' << v_s[i] << ' ' << rho[i] << '\n';
    }

    // Write the under_layers
    if (under_layers.size() > 0) {
        velocity_model << under_layers;
    }

    return velocity_model.str();
}


/**
 * Computes S- or P-wave first arrival
 * 
 * @author S.G. Solazzi
 * @author J. Cunha Teixeira
 * @param thk thickness of each layer [m]
 * @param vv velocity of each layer [m/s]
 * @param Xdata distance from the source to the receiver [m]
 * @param trig trigger time [s]
 * @return hodochrone table
*/
// 2021 - Solazzi - Implemented in Matlab
// MAR2024 - Cunha Teixeira - Traduced in C++
std::vector<double> firstArrivalSrc(const std::vector<double>& thk,
                                     const std::vector<double>& vv,
                                     const std::vector<double>& Xdata,
                                     const double& trig) {
    
    size_t rows = thk.size();
    size_t cols = Xdata.size();
    std::vector<double> Tr(rows*cols, 0.0); // arrival time table as a flatten array

    std::vector<double> Thod(cols, 0.0); // hodochrone table
    std::vector<double> Tz(rows, 0.0); // intercepts

    // Calculating intercepts
    for (size_t inl = 0; inl < rows; ++inl) {
        double Tzc = 0.0;
        for (size_t inltz = 0; inltz <= inl; ++inltz) {
            Tzc += 2 * thk[inltz] / vv[inltz] * sqrt(1 - vv[inltz] * vv[inltz] / (vv[inl] * vv[inl]));
        }
        Tz[inl] = Tzc;

        for (size_t ix = 0; ix < cols; ++ix) {
            Tr[inl + ix * rows] = Xdata[ix] / vv[inl] + Tz[inl];
        }
    }

    // Finding minimum arrival time for each Xdata point and adding trigger
    for (size_t ix = 0; ix < cols; ++ix) {
        double min_time = Tr[0 + ix * rows];
        for (size_t inl = 1; inl < rows; ++inl) {
            min_time = std::min(min_time, Tr[inl + ix * rows]);
        }
        Thod[ix] = min_time + trig;
    }

    return Thod;
}