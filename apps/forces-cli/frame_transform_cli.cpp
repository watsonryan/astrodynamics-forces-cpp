/**
 * @file frame_transform_cli.cpp
 * @brief Frame transform utility for external cross-validation (e.g., Astropy).
 * @author Watosn
 */

#include <cstdlib>
#include <string>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/core/transforms.hpp"

namespace {

double parse_double(const char* s) { return std::strtod(s, nullptr); }

void print_vec(const char* key, const astroforces::core::Vec3& v) {
  fmt::print("{}={:.17e},{:.17e},{:.17e}\n", key, v.x, v.y, v.z);
}

void print_usage() {
  spdlog::error(
      "usage:\n"
      "  frame_transform_cli gcrf_to_itrf <jd_utc> <jd_tt> <cip_x_rad> <cip_y_rad> <cip_s_rad> <xp_rad> <yp_rad> <dut1_s> "
      "<lod_s> <dX_rad> <dY_rad> <rx_gcrf_m> <ry_gcrf_m> <rz_gcrf_m> <vx_gcrf_mps> <vy_gcrf_mps> <vz_gcrf_mps> "
      "[cip_x_dot_rad_s] [cip_y_dot_rad_s] [cip_s_dot_rad_s] [xp_dot_rad_s] [yp_dot_rad_s] [dX_dot_rad_s] [dY_dot_rad_s]");
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    print_usage();
    return 1;
  }

  const std::string mode = argv[1];
  if (mode == "gcrf_to_itrf") {
    if (argc != 19 && argc != 26) {
      print_usage();
      return 3;
    }
    const double jd_utc = parse_double(argv[2]);
    const double jd_tt = parse_double(argv[3]);

    astroforces::core::CelestialIntermediatePole cip{};
    cip.x_rad = parse_double(argv[4]);
    cip.y_rad = parse_double(argv[5]);
    cip.s_rad = parse_double(argv[6]);

    astroforces::core::EarthOrientation eop{};
    eop.xp_rad = parse_double(argv[7]);
    eop.yp_rad = parse_double(argv[8]);
    eop.dut1_s = parse_double(argv[9]);
    eop.lod_s = parse_double(argv[10]);
    eop.dX_rad = parse_double(argv[11]);
    eop.dY_rad = parse_double(argv[12]);
    astroforces::core::CelestialIntermediatePoleRate cip_rate{};
    astroforces::core::EarthOrientationRate eop_rate{};
    const bool have_rates = (argc == 26);
    if (have_rates) {
      cip_rate.x_rad_s = parse_double(argv[19]);
      cip_rate.y_rad_s = parse_double(argv[20]);
      cip_rate.s_rad_s = parse_double(argv[21]);
      eop_rate.xp_rad_s = parse_double(argv[22]);
      eop_rate.yp_rad_s = parse_double(argv[23]);
      eop_rate.dX_rad_s = parse_double(argv[24]);
      eop_rate.dY_rad_s = parse_double(argv[25]);
    }

    const astroforces::core::Vec3 r_gcrf{
        parse_double(argv[13]), parse_double(argv[14]), parse_double(argv[15])};
    const astroforces::core::Vec3 v_gcrf{
        parse_double(argv[16]), parse_double(argv[17]), parse_double(argv[18])};

    astroforces::core::RotationWithDerivative rd{};
    if (have_rates) {
      rd = astroforces::core::gcrf_to_itrf_rotation_with_derivative_exact(jd_utc, jd_tt, cip, cip_rate, eop, eop_rate);
    } else {
      rd = astroforces::core::gcrf_to_itrf_rotation_with_derivative(jd_utc, jd_tt, cip, eop);
    }
    const auto r_itrf = astroforces::core::mat_vec(rd.r, r_gcrf);
    const auto v_itrf = astroforces::core::mat_vec(rd.r, v_gcrf) + astroforces::core::mat_vec(rd.dr, r_gcrf);
    const auto rt = astroforces::core::mat_transpose(rd.r);
    const auto r_gcrf_back = astroforces::core::mat_vec(rt, r_itrf);
    const auto v_gcrf_back = astroforces::core::mat_vec(rt, v_itrf - astroforces::core::mat_vec(rd.dr, r_gcrf_back));

    print_vec("r_itrf_m", r_itrf);
    print_vec("v_itrf_mps", v_itrf);
    print_vec("r_gcrf_back_m", r_gcrf_back);
    print_vec("v_gcrf_back_mps", v_gcrf_back);
    return 0;
  }

  print_usage();
  return 4;
}
