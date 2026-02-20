/**
 * @file test_third_body.cpp
 * @brief Third-body perturbation integration tests.
 * @author Watosn
 */

#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <iostream>

#include "dragcpp/atmo/types.hpp"
#include "dragcpp/forces/third_body.hpp"

namespace {

bool finite_vec(const dragcpp::atmo::Vec3& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

dragcpp::atmo::Vec3 add(const dragcpp::atmo::Vec3& a, const dragcpp::atmo::Vec3& b) {
  return dragcpp::atmo::Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
}

}  // namespace

int main() {
  namespace fs = std::filesystem;
  fs::path eph_path;
  if (const char* env = std::getenv("DRAGCPP_JPL_EPH_FILE")) {
    eph_path = env;
  } else {
    eph_path = fs::path(DRAGCPP_JPL_EPH_SOURCE_DIR) / "testdata" / "linux_p1550p2650.440";
  }
  if (!fs::exists(eph_path)) {
    std::cout << "third-body test skipped: ephemeris not found: " << eph_path << "\n";
    return 0;
  }

  dragcpp::atmo::StateVector state{};
  state.frame = dragcpp::atmo::Frame::ECI;
  state.epoch.utc_seconds = 1.0e9;
  state.position_m = dragcpp::atmo::Vec3{6778137.0, 0.0, 0.0};
  state.velocity_mps = dragcpp::atmo::Vec3{0.0, 7670.0, 0.0};

  const auto both = dragcpp::forces::ThirdBodyPerturbationModel::Create(
      {.ephemeris_file = eph_path, .use_sun = true, .use_moon = true, .name = "third_both"});
  const auto sun = dragcpp::forces::ThirdBodyPerturbationModel::Create(
      {.ephemeris_file = eph_path, .use_sun = true, .use_moon = false, .name = "third_sun"});
  const auto moon = dragcpp::forces::ThirdBodyPerturbationModel::Create(
      {.ephemeris_file = eph_path, .use_sun = false, .use_moon = true, .name = "third_moon"});

  const auto req = dragcpp::forces::PerturbationRequest{.state = state, .spacecraft = nullptr};
  const auto c_both = both->evaluate(req);
  const auto c_sun = sun->evaluate(req);
  const auto c_moon = moon->evaluate(req);
  if (c_both.status != dragcpp::atmo::Status::Ok || c_sun.status != dragcpp::atmo::Status::Ok ||
      c_moon.status != dragcpp::atmo::Status::Ok) {
    std::cerr << "third-body evaluation failed\n";
    return 1;
  }
  if (!finite_vec(c_both.acceleration_mps2) || !finite_vec(c_sun.acceleration_mps2) || !finite_vec(c_moon.acceleration_mps2)) {
    std::cerr << "third-body acceleration not finite\n";
    return 2;
  }
  if (!(dragcpp::atmo::norm(c_both.acceleration_mps2) > 0.0) || !(dragcpp::atmo::norm(c_sun.acceleration_mps2) > 0.0) ||
      !(dragcpp::atmo::norm(c_moon.acceleration_mps2) > 0.0)) {
    std::cerr << "third-body acceleration unexpectedly zero\n";
    return 3;
  }

  const auto sum = add(c_sun.acceleration_mps2, c_moon.acceleration_mps2);
  const auto diff = dragcpp::atmo::Vec3{c_both.acceleration_mps2.x - sum.x, c_both.acceleration_mps2.y - sum.y,
                                        c_both.acceleration_mps2.z - sum.z};
  if (!(dragcpp::atmo::norm(diff) <= 1e-16 * std::max(1.0, dragcpp::atmo::norm(c_both.acceleration_mps2)))) {
    std::cerr << "third-body superposition mismatch\n";
    return 4;
  }

  state.frame = dragcpp::atmo::Frame::ECEF;
  const auto c_bad_frame = both->evaluate(dragcpp::forces::PerturbationRequest{.state = state, .spacecraft = nullptr});
  if (c_bad_frame.status != dragcpp::atmo::Status::InvalidInput) {
    std::cerr << "expected invalid frame error\n";
    return 5;
  }

  return 0;
}
