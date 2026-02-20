/**
 * @file main.cpp
 * @brief astrodynamics-forces-cpp drag command-line entrypoint.
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
#include "dragcpp/weather/celestrak_csv_provider.hpp"
#include "dragcpp/weather/static_provider.hpp"

int main(int argc, char** argv) {
  if (argc < 8 || argc > 13) {
    std::cerr << "usage: drag_cli <x_m> <y_m> <z_m> <vx_mps> <vy_mps> <vz_mps> <epoch_utc_s> [model] [model_data] [wind] [wind_data] [weather_csv]\n";
    std::cerr << "models: basic | nrlmsis | dtm2020\n";
    std::cerr << "wind: zero | hwm14\n";
    std::cerr << "weather_csv: CelesTrak SW-Last5Years.csv (optional)\n";
    return 1;
  }

  astroforces::atmo::StateVector state{};
  state.position_m = astroforces::atmo::Vec3{std::atof(argv[1]), std::atof(argv[2]), std::atof(argv[3])};
  state.velocity_mps = astroforces::atmo::Vec3{std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6])};
  state.epoch.utc_seconds = std::atof(argv[7]);
  state.frame = astroforces::atmo::Frame::ECEF;

  const std::string weather_csv = (argc >= 13) ? argv[12] : "";
  const std::string model_name = (argc >= 9) ? argv[8] : "basic";
  const std::string model_data = (argc >= 10) ? argv[9] : "";
  const std::string wind_name = (argc >= 11) ? argv[10] : "zero";
const std::string wind_data = (argc >= 12) ? argv[11] : "";

  std::unique_ptr<astroforces::atmo::ISpaceWeatherProvider> weather{};
  if (!weather_csv.empty()) {
    weather = astroforces::weather::CelesTrakCsvSpaceWeatherProvider::Create(
        astroforces::weather::CelesTrakCsvSpaceWeatherProvider::Config{.csv_file = weather_csv});
  } else {
    const astroforces::atmo::WeatherIndices wx{.f107 = 150.0, .f107a = 150.0, .ap = 4.0, .kp = 2.0,
                                            .status = astroforces::atmo::Status::Ok};
    weather = std::make_unique<astroforces::weather::StaticSpaceWeatherProvider>(wx);
  }

  std::unique_ptr<astroforces::atmo::IAtmosphereModel> atmosphere{};
  if (model_name == "nrlmsis") {
    atmosphere = astroforces::adapters::Nrlmsis21AtmosphereAdapter::Create(
        astroforces::adapters::Nrlmsis21AtmosphereAdapter::Config{.parm_file = model_data});
  } else if (model_name == "dtm2020") {
    atmosphere = astroforces::adapters::Dtm2020AtmosphereAdapter::Create(
        astroforces::adapters::Dtm2020AtmosphereAdapter::Config{.coeff_file = model_data});
  } else {
    atmosphere = std::make_unique<astroforces::models::ExponentialAtmosphereModel>(1.225, 0.0, 7000.0, 1000.0);
  }
  std::unique_ptr<astroforces::atmo::IWindModel> wind{};
  if (wind_name == "hwm14") {
    wind = astroforces::adapters::Hwm14WindAdapter::Create(
        astroforces::adapters::Hwm14WindAdapter::Config{.data_dir = wind_data});
  } else {
    wind = std::make_unique<astroforces::models::ZeroWindModel>();
  }

  astroforces::sc::SpacecraftProperties sc{.mass_kg = 1200.0,
                                       .reference_area_m2 = 12.0,
                                       .cd = 2.2,
                                       .use_surface_model = false,
                                       .surfaces = {}};

  astroforces::drag::DragAccelerationModel model(*weather, *atmosphere, *wind);
  const auto result = model.evaluate(state, sc);

  if (result.status != astroforces::atmo::Status::Ok) {
    std::cerr << "drag eval failed\n";
    return 2;
  }

  std::cout << "ax=" << result.acceleration_mps2.x << " ay=" << result.acceleration_mps2.y
            << " az=" << result.acceleration_mps2.z << "\n";
  std::cout << "rho=" << result.density_kg_m3 << " temp_k=" << result.temperature_k << " vrel_mps=" << result.relative_speed_mps
            << " q_pa=" << result.dynamic_pressure_pa << " area=" << result.area_m2 << " cd_eff=" << result.cd << "\n";
  std::cout << "wx_source=" << static_cast<int>(result.weather.source) << " wx_interp=" << (result.weather.interpolated ? 1 : 0)
            << " wx_extrap=" << (result.weather.extrapolated ? 1 : 0) << " f107=" << result.weather.f107
            << " f107a=" << result.weather.f107a << " ap=" << result.weather.ap << " kp=" << result.weather.kp
            << " ap3h=" << result.weather.ap_3h_current << " kp3h=" << result.weather.kp_3h_current << "\n";
  return 0;
}
