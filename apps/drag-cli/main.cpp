/**
 * @file main.cpp
 * @brief drag-cpp command-line entrypoint.
 * @author Watosn
 */

#include <cstdlib>
#include <iostream>

#include "dragcpp/drag/drag_model.hpp"
#include "dragcpp/models/exponential_atmosphere.hpp"
#include "dragcpp/sc/spacecraft.hpp"
#include "dragcpp/weather/static_provider.hpp"

int main(int argc, char** argv) {
  if (argc != 8) {
    std::cerr << "usage: drag_cli <x_m> <y_m> <z_m> <vx_mps> <vy_mps> <vz_mps> <epoch_utc_s>\n";
    return 1;
  }

  dragcpp::atmo::StateVector state{};
  state.position_m = dragcpp::atmo::Vec3{std::atof(argv[1]), std::atof(argv[2]), std::atof(argv[3])};
  state.velocity_mps = dragcpp::atmo::Vec3{std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6])};
  state.epoch.utc_seconds = std::atof(argv[7]);
  state.frame = dragcpp::atmo::Frame::ECI;

  const dragcpp::atmo::WeatherIndices wx{.f107 = 150.0, .f107a = 150.0, .ap = 4.0, .kp = 2.0,
                                          .status = dragcpp::atmo::Status::Ok};
  dragcpp::weather::StaticSpaceWeatherProvider weather(wx);
  dragcpp::models::ExponentialAtmosphereModel atmosphere(1.225, 0.0, 7000.0, 1000.0);
  dragcpp::models::ZeroWindModel wind;

  dragcpp::sc::SpacecraftProperties sc{.mass_kg = 1200.0,
                                       .reference_area_m2 = 12.0,
                                       .cd = 2.2,
                                       .use_surface_model = false,
                                       .surfaces = {}};

  dragcpp::drag::DragAccelerationModel model(weather, atmosphere, wind);
  const auto result = model.evaluate(state, sc);

  if (result.status != dragcpp::atmo::Status::Ok) {
    std::cerr << "drag eval failed\n";
    return 2;
  }

  std::cout << "ax=" << result.acceleration_mps2.x << " ay=" << result.acceleration_mps2.y
            << " az=" << result.acceleration_mps2.z << "\n";
  std::cout << "rho=" << result.density_kg_m3 << " area=" << result.area_m2 << " cd=" << result.cd << "\n";
  return 0;
}
