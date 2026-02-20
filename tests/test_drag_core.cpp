/**
 * @file test_drag_core.cpp
 * @brief Core drag equation tests.
 * @author Watosn
 */

#include <cmath>
#include <iostream>

#include "dragcpp/drag/drag_model.hpp"
#include "dragcpp/models/exponential_atmosphere.hpp"
#include "dragcpp/sc/spacecraft.hpp"
#include "dragcpp/weather/static_provider.hpp"

namespace {

bool approx(double a, double b, double rel) {
  const double d = std::abs(a - b);
  const double n = std::max(std::abs(b), 1e-30);
  return d / n <= rel;
}

}  // namespace

int main() {
  using namespace dragcpp;
  const atmo::WeatherIndices wx{.f107 = 150.0, .f107a = 150.0, .ap = 4.0, .kp = 2.0, .status = atmo::Status::Ok};
  weather::StaticSpaceWeatherProvider weather(wx);
  models::ExponentialAtmosphereModel atmosphere(1.225, 0.0, 7000.0, 1000.0);
  models::ZeroWindModel wind;
  drag::DragAccelerationModel model(weather, atmosphere, wind);

  atmo::StateVector state{};
  state.frame = atmo::Frame::ECI;
  state.position_m = atmo::Vec3{6378137.0, 0.0, 0.0};
  state.velocity_mps = atmo::Vec3{7500.0, 0.0, 0.0};

  sc::SpacecraftProperties sc{.mass_kg = 1000.0, .reference_area_m2 = 10.0, .cd = 2.0, .use_surface_model = false};

  const auto out = model.evaluate(state, sc);
  if (out.status != atmo::Status::Ok) {
    std::cerr << "status failed\n";
    return 1;
  }

  const double expected_ax = -0.5 * 1.225 * 2.0 * 10.0 / 1000.0 * 7500.0 * 7500.0;
  if (!approx(out.acceleration_mps2.x, expected_ax, 1e-12)) {
    std::cerr << "ax mismatch\n";
    return 2;
  }
  if (!approx(out.acceleration_mps2.y, 0.0, 1e-12) || !approx(out.acceleration_mps2.z, 0.0, 1e-12)) {
    std::cerr << "vector mismatch\n";
    return 3;
  }

  return 0;
}
