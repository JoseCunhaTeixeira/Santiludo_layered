#include <cmath>
#include <stdexcept>
#include <vector>
#include <string>
#include <sstream>
#include <omp.h>
#include <algorithm>



std::string writeVelocityModel_src(const std::vector<double>& thk,
                        const std::vector<double>& vp,
                        const std::vector<double>& vs,
                        const std::vector<double>& rho,
                        const std::string under_layers,
                        const int n_under_layers) {
    // ==========================================
    // Returns a string of the velocity model in a format expected by GPDC
    //
    // Parameters:
    // vector<double> thk: thickness of each layer [m]
    // vector<double> vp: P-wave velocity of each layer [m/s]
    // vector<double> vs: S-wave velocity of each layer [m/s]
    // vector<double> rho: density of each layer [kg/m^3]
    // string under_layers: Layers to put under the studied soil column on the velocity model in GPDC format
    // int n_under_layers: Number of layers in the under_layers
    //
    // Outputs:
    // string velocity_model: velocity model in GPDC format
    //
    // Programmers:
    // Firstly implemented in Matlab by S. Pasquet in Solazzi et al. (2021)
    // Traduced in C++ by J. Cunha Teixeira in 2024/03
    // ==========================================

    int nl = thk.size(); // Number of layers
    std::stringstream velocity_model; // String to store the velocity model
    velocity_model << nl+n_under_layers << '\n'; // Number of layers in the velocity model

    // Write thickness and velocity in dinver format
    for (int i = 0; i < nl; ++i) {
        if (i == nl - 1 && under_layers.size() == 0) {
            velocity_model << "0 ";
        } else {
            velocity_model << thk[i] << ' ';
        }
        velocity_model << vp[i] << ' ' << vs[i] << ' ' << rho[i] << '\n';
    }

    // Write the under_layers
    if (under_layers.size() > 0) {
        velocity_model << under_layers;
    }

    return velocity_model.str();
}



// WITH PARALLELIZATION - Must be compiled with '-fopenmp' in setup.py and needs #include <omp.h>
std::vector<double> firstArrival_src(const std::vector<double>& thk,
                                     const std::vector<double>& vv,
                                     const std::vector<double>& Xdata,
                                     double trig) {
    // ==========================================
    // Computes S- or P-wave first arrivals
    //
    // Parameters:
    // vector<double> thk: thickness of each layer [m]
    // vector<double> vv: velocity of each layer [m/s]
    // vector<double> Xdata: distance from the source to the receiver [m]
    // double trig: trigger time [s]
    //
    // Outputs:
    // vector<double> Thod: hodochrone table
    //
    // Programmers:
    // Firstly implemented in Matlab by L. Bodet in Solazzi et al. (2021)
    // Traduced in C++ by J. Cunha Teixeira in 2024/03
    // ==========================================
    
    std::vector<double> Thod(Xdata.size(), 0.0); // Hodochrone table
    std::vector<double> Tz(thk.size(), 0.0); // Intercepts

    // Calculating intercepts in parallel
#pragma omp parallel for
    for (int inl = 0; inl < static_cast<int>(thk.size()); ++inl) {
        double Tzc = 0.0;
        for (size_t inltz = 0; inltz <= static_cast<size_t>(inl); ++inltz) {
            Tzc += 2 * thk[inltz] / vv[inltz] * sqrt(1 - vv[inltz] * vv[inltz] / (vv[inl] * vv[inl]));
        }
        Tz[inl] = Tzc;
    }

    // Finding minimum arrival time for each Xdata point and adding trigger
#pragma omp parallel for
    for (size_t ix = 0; ix < Xdata.size(); ++ix) {
        double min_time = Xdata[ix] / vv[0] + Tz[0];
        for (size_t inl = 1; inl < thk.size(); ++inl) {
            double arrival_time = Xdata[ix] / vv[inl] + Tz[inl];
            min_time = (arrival_time < min_time) ? arrival_time : min_time;
        }
        Thod[ix] = min_time + trig;
    }

    return Thod;
}



// WITHOUT PARALLELIZATION - Do not need to be compiled with '-fopenmp' in setup.py and #include <omp.h>
// std::vector<double> firstArrival_src(const std::vector<double>& thk,
//                                      const std::vector<double>& vv,
//                                      const std::vector<double>& Xdata,
//                                      double trig) {
//     // ==========================================
//     // Computes S- or P-wave first arrivals
//     //
//     // Parameters:
//     // vector<double> thk: thickness of each layer [m]
//     // vector<double> vv: velocity of each layer [m/s]
//     // vector<double> Xdata: distance from the source to the receiver [m]
//     // double trig: trigger time [s]
//     //
//     // Outputs:
//     // vector<double> Thod: hodochrone table
//     //
//     // Programmers:
//     // Firstly implemented in Matlab by L. Bodet in Solazzi et al. (2021)
//     // Traduced in C++ by J. Cunha Teixeira in 2024/03
//     // ==========================================

//     std::vector<std::vector<double>> Tr(thk.size(), std::vector<double>(Xdata.size(), 0.0)); // arrival time table
//     std::vector<double> Thod(Xdata.size(), 0.0); // hodochrone table
//     std::vector<double> Tz(thk.size(), 0.0); // intercepts

//     // Calculating intercepts
//     for (size_t inl = 0; inl < thk.size(); ++inl) {
//         double Tzc = 0.0;
//         for (size_t inltz = 0; inltz <= inl; ++inltz) {
//             Tzc += 2 * thk[inltz] / vv[inltz] * sqrt(1 - vv[inltz] * vv[inltz] / (vv[inl] * vv[inl]));
//         }
//         Tz[inl] = Tzc;

//         for (size_t ix = 0; ix < Xdata.size(); ++ix) {
//             Tr[inl][ix] = Xdata[ix] / vv[inl] + Tz[inl];
//         }
//     }

//     // Finding minimum arrival time for each Xdata point and adding trigger
//     for (size_t ix = 0; ix < Xdata.size(); ++ix) {
//         double min_time = Tr[0][ix];
//         for (size_t inl = 1; inl < thk.size(); ++inl) {
//             min_time = std::min(min_time, Tr[inl][ix]);
//         }
//         Thod[ix] = min_time + trig;
//     }

//     return Thod;
// }