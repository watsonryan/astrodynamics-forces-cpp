/**
 * @file pole_tide.cpp
 * @brief Pole tide correction helpers for gravity SH models.
 * @author Watosn
 */

#include "astroforces/forces/gravity/tides/pole_tide.hpp"

#include "astroforces/atmo/constants.hpp"

namespace astroforces::forces::tides {
namespace {

void secular_pole_mas(const double mjd_tt, double& xpv_mas, double& ypv_mas) {
  const double t_years = (mjd_tt - 51544.5) / 365.25;
  xpv_mas = 55.0 + 1.677 * t_years;
  ypv_mas = 320.5 + 3.460 * t_years;
}

}  // namespace

void add_pole_solid_tide_delta(const double mjd_tt,
                               const double xp_rad,
                               const double yp_rad,
                               Eigen::MatrixXd& dC,
                               Eigen::MatrixXd& dS) {
  if (dC.rows() <= 2 || dC.cols() <= 1 || dS.rows() <= 2 || dS.cols() <= 1) {
    return;
  }

  double xpv_mas = 0.0;
  double ypv_mas = 0.0;
  secular_pole_mas(mjd_tt, xpv_mas, ypv_mas);

  const double m1 = +(xp_rad / astroforces::core::constants::kArcsecToRad - xpv_mas / 1000.0);
  const double m2 = -(yp_rad / astroforces::core::constants::kArcsecToRad - ypv_mas / 1000.0);

  dC(2, 1) += -1.333e-9 * (m1 + 0.0115 * m2);
  dS(2, 1) += -1.333e-9 * (m2 - 0.0115 * m1);
}

void add_pole_ocean_tide_delta(const double mjd_tt,
                               const double xp_rad,
                               const double yp_rad,
                               Eigen::MatrixXd& dC,
                               Eigen::MatrixXd& dS) {
  if (dC.rows() <= 2 || dC.cols() <= 1 || dS.rows() <= 2 || dS.cols() <= 1) {
    return;
  }

  double xpv_mas = 0.0;
  double ypv_mas = 0.0;
  secular_pole_mas(mjd_tt, xpv_mas, ypv_mas);

  const double m1 = +(xp_rad / astroforces::core::constants::kArcsecToRad - xpv_mas / 1000.0);
  const double m2 = -(yp_rad / astroforces::core::constants::kArcsecToRad - ypv_mas / 1000.0);

  dC(2, 1) += -2.1778e-10 * (m1 - 0.01724 * m2);
  dS(2, 1) += -1.7232e-10 * (m2 - 0.03365 * m1);
}

}  // namespace astroforces::forces::tides
