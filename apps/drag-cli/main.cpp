/**
 * @file main.cpp
 * @brief drag-cpp command-line entrypoint.
 * @author Watosn
 */

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

#include "dragcpp/adapters/dtm2020_adapter.hpp"
#include "dragcpp/adapters/hwm14_adapter.hpp"
#include "dragcpp/adapters/nrlmsis21_adapter.hpp"
#include "dragcpp/drag/drag_model.hpp"
#include "dragcpp/models/exponential_atmosphere.hpp"
#include "dragcpp/sc/spacecraft.hpp"
#include "dragcpp/weather/static_provider.hpp"

int main(int argc, char** argv) {
  if (argc < 8 || argc > 12) {
    std::cerr << "usage: drag_cli <x_m> <y_m> <z_m> <vx_mps> <vy_mps> <vz_mps> <epoch_utc_s> [model] [model_data] [wind] [wind_data]\n";
    std::cerr << "models: basic | nrlmsis | dtm2020\n";
    std::cerr << "wind: zero | hwm14\n";
    return 1;
  }

  dragcpp::atmo::StateVector state{};
  state.position_m = dragcpp::atmo::Vec3{std::atof(argv[1]), std::atof(argv[2]), std::atof(argv[3])};
  state.velocity_mps = dragcpp::atmo::Vec3{std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6])};
  state.epoch.utc_seconds = std::atof(argv[7]);
  state.frame = dragcpp::atmo::Frame::ECEF;

  const dragcpp::atmo::WeatherIndices wx{.f107 = 150.0, .f107a = 150.0, .ap = 4.0, .kp = 2.0,
                                          .status = dragcpp::atmo::Status::Ok};
  dragcpp::weather::StaticSpaceWeatherProvider weather(wx);
  const std::string model_name = (argc >= 9) ? argv[8] : "basic";
  const std::string model_data = (argc >= 10) ? argv[9] : "";
  const std::string wind_name = (argc >= 11) ? argv[10] : "zero";
  const std::string wind_data = (argc >= 12) ? argv[11] : "";

  std::unique_ptr<dragcpp::atmo::IAtmosphereModel> atmosphere{};
  if (model_name == "nrlmsis") {
    atmosphere = dragcpp::adapters::Nrlmsis21AtmosphereAdapter::Create(
        dragcpp::adapters::Nrlmsis21AtmosphereAdapter::Config{.parm_file = model_data});
  } else if (model_name == "dtm2020") {
    atmosphere = dragcpp::adapters::Dtm2020AtmosphereAdapter::Create(
        dragcpp::adapters::Dtm2020AtmosphereAdapter::Config{.coeff_file = model_data});
  } else {
    atmosphere = std::make_unique<dragcpp::models::ExponentialAtmosphereModel>(1.225, 0.0, 7000.0, 1000.0);
  }
  std::unique_ptr<dragcpp::atmo::IWindModel> wind{};
  if (wind_name == "hwm14") {
    wind = dragcpp::adapters::Hwm14WindAdapter::Create(
        dragcpp::adapters::Hwm14WindAdapter::Config{.data_dir = wind_data});
  } else {
    wind = std::make_unique<dragcpp::models::ZeroWindModel>();
  }

  dragcpp::sc::SpacecraftProperties sc{.mass_kg = 1200.0,
                                       .reference_area_m2 = 12.0,
                                       .cd = 2.2,
                                       .use_surface_model = false,
                                       .surfaces = {}};

  dragcpp::drag::DragAccelerationModel model(weather, *atmosphere, *wind);
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
