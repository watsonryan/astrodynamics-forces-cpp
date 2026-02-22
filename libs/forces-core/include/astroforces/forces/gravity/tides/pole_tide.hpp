/**
 * @file pole_tide.hpp
 * @brief Pole tide correction helpers for gravity SH models.
 * @author Watosn
 */
#pragma once

#include <Eigen/Dense>

namespace astroforces::forces::tides {

void add_pole_solid_tide_delta(double mjd_tt,
                               double xp_rad,
                               double yp_rad,
                               Eigen::MatrixXd& dC,
                               Eigen::MatrixXd& dS);

void add_pole_ocean_tide_delta(double mjd_tt,
                               double xp_rad,
                               double yp_rad,
                               Eigen::MatrixXd& dC,
                               Eigen::MatrixXd& dS);

}  // namespace astroforces::forces::tides
