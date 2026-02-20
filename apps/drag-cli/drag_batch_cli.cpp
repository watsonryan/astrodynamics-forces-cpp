/**
 * @file drag_batch_cli.cpp
 * @brief Batch drag profile runner with CSV/JSON outputs.
 * @author Watosn
 */

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <ctime>
#include <vector>

#include <spdlog/spdlog.h>

#include "dragcpp/adapters/dtm2020_adapter.hpp"
#include "dragcpp/adapters/hwm14_adapter.hpp"
#include "dragcpp/adapters/nrlmsis21_adapter.hpp"
#include "dragcpp/drag/drag_model.hpp"
#include "dragcpp/models/exponential_atmosphere.hpp"
#include "dragcpp/sc/spacecraft.hpp"
#include "dragcpp/weather/celestrak_csv_provider.hpp"
#include "dragcpp/weather/static_provider.hpp"

namespace {

const char* status_to_string(astroforces::atmo::Status s) {
  switch (s) {
    case astroforces::atmo::Status::Ok:
      return "ok";
    case astroforces::atmo::Status::InvalidInput:
      return "invalid_input";
    case astroforces::atmo::Status::NotImplemented:
      return "not_implemented";
    case astroforces::atmo::Status::DataUnavailable:
      return "data_unavailable";
    case astroforces::atmo::Status::NumericalError:
      return "numerical_error";
    default:
      return "unknown";
  }
}

const char* weather_source_to_string(astroforces::atmo::WeatherSource s) {
  switch (s) {
    case astroforces::atmo::WeatherSource::StaticProvider:
      return "static";
    case astroforces::atmo::WeatherSource::CelesTrakLast5YearsCsv:
      return "celestrak_last5y";
    default:
      return "unknown";
  }
}

struct SampleRow {
  double epoch_utc_s{};
  astroforces::atmo::Vec3 position_m{};
  astroforces::atmo::Vec3 velocity_mps{};
};

bool parse_sample_row(const std::string& line, SampleRow& out) {
  std::stringstream ss(line);
  std::string tok;
  std::vector<double> values;
  while (std::getline(ss, tok, ',')) {
    if (tok.empty()) {
      return false;
    }
    char* end = nullptr;
    const double v = std::strtod(tok.c_str(), &end);
    if (end == tok.c_str() || *end != '\0') {
      return false;
    }
    values.push_back(v);
  }
  if (values.size() != 7U) {
    return false;
  }
  out.epoch_utc_s = values[0];
  out.position_m = astroforces::atmo::Vec3{values[1], values[2], values[3]};
  out.velocity_mps = astroforces::atmo::Vec3{values[4], values[5], values[6]};
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 4 || argc > 9) {
    spdlog::error("usage: drag_batch_cli <input_csv> <output_file> <format:csv|json> [model] [model_data] [wind] [wind_data] [weather_csv]");
    spdlog::error("input row: epoch_utc_s,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps");
    return 1;
  }

  const std::filesystem::path input_path = argv[1];
  const std::filesystem::path output_path = argv[2];
  const std::string format = argv[3];
  const std::string model_name = (argc >= 5) ? argv[4] : "basic";
  const std::string model_data = (argc >= 6) ? argv[5] : "";
  const std::string wind_name = (argc >= 7) ? argv[6] : "zero";
  const std::string wind_data = (argc >= 8) ? argv[7] : "";
  const std::string weather_csv = (argc >= 9) ? argv[8] : "";
  if (format != "csv" && format != "json") {
    spdlog::error("format must be csv or json");
    return 4;
  }

  std::ifstream in(input_path);
  if (!in) {
    spdlog::error("failed to open input csv: {}", input_path.string());
    return 2;
  }
  std::ofstream out(output_path);
  if (!out) {
    spdlog::error("failed to open output file: {}", output_path.string());
    return 3;
  }

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

  astroforces::sc::SpacecraftProperties sc{
      .mass_kg = 1200.0, .reference_area_m2 = 12.0, .cd = 2.2, .use_surface_model = false, .surfaces = {}};
  astroforces::drag::DragAccelerationModel drag(*weather, *atmosphere, *wind);

  std::time_t now = std::time(nullptr);
  if (format == "csv") {
    out << "#record_type=metadata,schema=drag_batch_v1,project=astrodynamics-forces-cpp,generated_unix_utc=" << now
        << ",model=" << model_name << ",model_data=" << model_data << ",wind=" << wind_name << ",wind_data=" << wind_data
        << ",weather_csv=" << weather_csv << "\n";
    out << "epoch_utc_s,ax_mps2,ay_mps2,az_mps2,rho_kg_m3,temp_k,vrel_mps,q_pa,area_m2,cd_eff,"
           "wx_source,wx_interp,wx_extrap,f107,f107a,ap_daily,kp_daily,ap_3h,kp_3h,status\n";
  } else {
    out << "{\"record_type\":\"metadata\",\"schema\":\"drag_batch_v1\",\"project\":\"astrodynamics-forces-cpp\","
           "\"generated_unix_utc\":" << now << ",\"model\":\"" << model_name << "\",\"model_data\":\"" << model_data
        << "\",\"wind\":\"" << wind_name << "\",\"wind_data\":\"" << wind_data << "\",\"weather_csv\":\"" << weather_csv
        << "\"}\n";
  }

  std::string line;
  std::size_t line_no = 0;
  while (std::getline(in, line)) {
    ++line_no;
    if (line.empty()) {
      continue;
    }
    SampleRow row{};
    if (!parse_sample_row(line, row)) {
      if (line_no == 1 && line.find("epoch_utc_s") != std::string::npos) {
        continue;
      }
      spdlog::warn("skipping malformed row {}", line_no);
      continue;
    }

    astroforces::atmo::StateVector state{};
    state.epoch.utc_seconds = row.epoch_utc_s;
    state.position_m = row.position_m;
    state.velocity_mps = row.velocity_mps;
    state.frame = astroforces::atmo::Frame::ECEF;

    const auto r = drag.evaluate(state, sc);
    if (format == "json") {
      out << "{\"record_type\":\"sample\",\"schema\":\"drag_batch_v1\",\"epoch_utc_s\":" << row.epoch_utc_s
          << ",\"ax_mps2\":" << r.acceleration_mps2.x
          << ",\"ay_mps2\":" << r.acceleration_mps2.y << ",\"az_mps2\":" << r.acceleration_mps2.z
          << ",\"rho_kg_m3\":" << r.density_kg_m3 << ",\"temp_k\":" << r.temperature_k
          << ",\"vrel_mps\":" << r.relative_speed_mps << ",\"q_pa\":" << r.dynamic_pressure_pa
          << ",\"area_m2\":" << r.area_m2 << ",\"cd_eff\":" << r.cd << ",\"wx_source\":\""
          << weather_source_to_string(r.weather.source) << "\",\"wx_interp\":" << (r.weather.interpolated ? "true" : "false")
          << ",\"wx_extrap\":" << (r.weather.extrapolated ? "true" : "false") << ",\"f107\":" << r.weather.f107
          << ",\"f107a\":" << r.weather.f107a << ",\"ap_daily\":" << r.weather.ap << ",\"kp_daily\":" << r.weather.kp
          << ",\"ap_3h\":" << r.weather.ap_3h_current << ",\"kp_3h\":" << r.weather.kp_3h_current << ",\"status\":\""
          << status_to_string(r.status) << "\"}\n";
    } else {
      out << row.epoch_utc_s << "," << r.acceleration_mps2.x << "," << r.acceleration_mps2.y << ","
          << r.acceleration_mps2.z << "," << r.density_kg_m3 << "," << r.temperature_k << ","
          << r.relative_speed_mps << "," << r.dynamic_pressure_pa << "," << r.area_m2 << "," << r.cd << ","
          << weather_source_to_string(r.weather.source) << "," << (r.weather.interpolated ? 1 : 0) << ","
          << (r.weather.extrapolated ? 1 : 0) << "," << r.weather.f107 << "," << r.weather.f107a << ","
          << r.weather.ap << "," << r.weather.kp << "," << r.weather.ap_3h_current << ","
          << r.weather.kp_3h_current << "," << status_to_string(r.status) << "\n";
    }
  }

  return 0;
}
