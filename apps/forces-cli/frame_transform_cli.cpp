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
      "  frame_transform_cli simple <utc_s> <rx_eci_m> <ry_eci_m> <rz_eci_m> <vx_eci_mps> <vy_eci_mps> <vz_eci_mps>\n"
      "  frame_transform_cli gcrf_to_itrf <jd_utc> <jd_tt> <cip_x_rad> <cip_y_rad> <cip_s_rad> <xp_rad> <yp_rad> <dut1_s> "
      "<lod_s> <dX_rad> <dY_rad> <rx_gcrf_m> <ry_gcrf_m> <rz_gcrf_m> <vx_gcrf_mps> <vy_gcrf_mps> <vz_gcrf_mps>");
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    print_usage();
    return 1;
  }

  const std::string mode = argv[1];
  if (mode == "simple") {
    if (argc != 9) {
      print_usage();
      return 2;
    }
    const double utc_s = parse_double(argv[2]);
    const astroforces::core::Vec3 r_eci{
        parse_double(argv[3]), parse_double(argv[4]), parse_double(argv[5])};
    const astroforces::core::Vec3 v_eci{
        parse_double(argv[6]), parse_double(argv[7]), parse_double(argv[8])};

    const auto r_ecef = astroforces::core::eci_to_ecef_position(r_eci, utc_s);
    const auto v_ecef = astroforces::core::eci_to_ecef_velocity(r_eci, v_eci, utc_s);
    const auto r_eci_back = astroforces::core::ecef_to_eci_position(r_ecef, utc_s);
    const auto v_eci_back = astroforces::core::ecef_to_eci_velocity(r_ecef, v_ecef, utc_s);

    print_vec("r_ecef_m", r_ecef);
    print_vec("v_ecef_mps", v_ecef);
    print_vec("r_eci_back_m", r_eci_back);
    print_vec("v_eci_back_mps", v_eci_back);
    return 0;
  }

  if (mode == "gcrf_to_itrf") {
    if (argc != 19) {
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

    const astroforces::core::Vec3 r_gcrf{
        parse_double(argv[13]), parse_double(argv[14]), parse_double(argv[15])};
    const astroforces::core::Vec3 v_gcrf{
        parse_double(argv[16]), parse_double(argv[17]), parse_double(argv[18])};

    const auto r_itrf = astroforces::core::gcrf_to_itrf_position(r_gcrf, jd_utc, jd_tt, cip, eop);
    const auto v_itrf = astroforces::core::gcrf_to_itrf_velocity(r_gcrf, v_gcrf, jd_utc, jd_tt, cip, eop);
    const auto r_gcrf_back = astroforces::core::itrf_to_gcrf_position(r_itrf, jd_utc, jd_tt, cip, eop);
    const auto v_gcrf_back = astroforces::core::itrf_to_gcrf_velocity(r_itrf, v_itrf, jd_utc, jd_tt, cip, eop);

    print_vec("r_itrf_m", r_itrf);
    print_vec("v_itrf_mps", v_itrf);
    print_vec("r_gcrf_back_m", r_gcrf_back);
    print_vec("v_gcrf_back_mps", v_gcrf_back);
    return 0;
  }

  print_usage();
  return 4;
}
