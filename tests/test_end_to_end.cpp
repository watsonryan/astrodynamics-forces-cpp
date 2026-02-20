/**
 * @file test_end_to_end.cpp
 * @brief End-to-end drag pipeline smoke test.
 * @author Watosn
 */

#include <iostream>

#include "dragcpp/drag/drag_model.hpp"
#include "dragcpp/models/exponential_atmosphere.hpp"
#include "dragcpp/sc/spacecraft.hpp"
#include "dragcpp/weather/static_provider.hpp"

int main() {
  using namespace astroforces;

  const atmo::WeatherIndices wx{.f107 = 120.0, .f107a = 130.0, .ap = 8.0, .kp = 3.0, .status = atmo::Status::Ok};
  weather::StaticSpaceWeatherProvider weather(wx);
  models::ExponentialAtmosphereModel atmosphere(3.0e-11, 400e3, 65e3, 900.0);
  models::ZeroWindModel wind;
  drag::DragAccelerationModel model(weather, atmosphere, wind);

  atmo::StateVector state{};
  state.epoch.utc_seconds = 1.0e9;
  state.frame = atmo::Frame::ECI;
  state.position_m = atmo::Vec3{6778137.0, 0.0, 0.0};
  state.velocity_mps = atmo::Vec3{0.0, 7670.0, 0.0};

  sc::SpacecraftProperties sc{.mass_kg = 600.0, .reference_area_m2 = 4.0, .cd = 2.25, .use_surface_model = false};

  const auto out = model.evaluate(state, sc);
  if (out.status != atmo::Status::Ok) {
    std::cerr << "pipeline failed\n";
    return 1;
  }
  if (!(out.density_kg_m3 > 0.0)) {
    std::cerr << "density invalid\n";
    return 2;
  }
  return 0;
}
