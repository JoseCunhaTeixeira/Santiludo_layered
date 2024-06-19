#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>
#include <string>

struct vanGenResult {
    std::vector<double> h, s_w, s_we;
};

struct soilType {
    double wsand, wclay, wsilt, phi, alpha, nvg, theta, s_wr;
};

soilType selectSoilTypeSrc(const std::string_view& soil_type);

vanGenResult vanGenSrc(const std::vector<double>& z,
                         const double& wt,
                         const std::vector<std::string>& soil_types,
                         const std::vector<double>& thicknesses);