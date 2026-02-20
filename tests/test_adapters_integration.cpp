/**
 * @file test_adapters_integration.cpp
 * @brief Adapter-backed integration tests across NRLMSIS, DTM2020, and HWM14.
 * @author Watosn
 */

#include <cmath>
#include <filesystem>
#include <iostream>

#include "dragcpp/adapters/dtm2020_adapter.hpp"
#include "dragcpp/adapters/hwm14_adapter.hpp"
#include "dragcpp/adapters/nrlmsis21_adapter.hpp"
#include "dragcpp/drag/drag_model.hpp"
#include "dragcpp/sc/spacecraft.hpp"
#include "dragcpp/weather/static_provider.hpp"

namespace {

bool finite(const dragcpp::atmo::Vec3& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

}  // namespace

int main() {
  using namespace dragcpp;
  namespace fs = std::filesystem;

  const auto nrl_parm = fs::path(DRAGCPP_NRLMSIS21_SOURCE_DIR) / "data" / "msis21.parm";
  const auto dtm_coeff = fs::path(DRAGCPP_DTM2020_SOURCE_DIR) / "testdata" / "operational_regression_coeff.dat";
  const auto hwm_data_dir = fs::path(DRAGCPP_HWM14_SOURCE_DIR) / "testdata";

  if (!fs::exists(nrl_parm) || !fs::exists(dtm_coeff) || !fs::exists(hwm_data_dir)) {
    std::cerr << "adapter test data paths missing\n";
    return 10;
  }

  const atmo::WeatherIndices wx{.f107 = 150.0, .f107a = 150.0, .ap = 4.0, .kp = 2.0, .status = atmo::Status::Ok};

  atmo::StateVector state{};
  state.epoch.utc_seconds = 1.0e9;
  state.frame = atmo::Frame::ECEF;
  state.position_m = atmo::Vec3{6778137.0, 0.0, 0.0};
  state.velocity_mps = atmo::Vec3{0.0, 7670.0, 0.0};

  const auto nrl = adapters::Nrlmsis21AtmosphereAdapter::Create({.parm_file = nrl_parm});
  const auto nrl_out = nrl->evaluate(state, wx);
  if (nrl_out.status != atmo::Status::Ok || !(nrl_out.density_kg_m3 > 0.0) || !std::isfinite(nrl_out.temperature_k)) {
    std::cerr << "nrlmsis adapter evaluation failed\n";
    return 1;
  }

  const auto dtm = adapters::Dtm2020AtmosphereAdapter::Create({.coeff_file = dtm_coeff});
  const auto dtm_out = dtm->evaluate(state, wx);
  if (dtm_out.status != atmo::Status::Ok || !(dtm_out.density_kg_m3 > 0.0) || !std::isfinite(dtm_out.temperature_k)) {
    std::cerr << "dtm2020 adapter evaluation failed\n";
    return 2;
  }

  const auto hwm = adapters::Hwm14WindAdapter::Create({.data_dir = hwm_data_dir});
  const auto hwm_out = hwm->evaluate(state, wx);
  if (hwm_out.status != atmo::Status::Ok || hwm_out.frame != atmo::Frame::ECEF || !finite(hwm_out.velocity_mps)) {
    std::cerr << "hwm14 adapter evaluation failed\n";
    return 3;
  }

  weather::StaticSpaceWeatherProvider weather(wx);
  sc::SpacecraftProperties sc{.mass_kg = 600.0, .reference_area_m2 = 4.0, .cd = 2.25, .use_surface_model = false};
  drag::DragAccelerationModel drag_model(weather, *nrl, *hwm);
  const auto drag_out = drag_model.evaluate(state, sc);
  if (drag_out.status != atmo::Status::Ok || !(drag_out.density_kg_m3 > 0.0) || !finite(drag_out.acceleration_mps2)) {
    std::cerr << "adapter-backed drag pipeline failed\n";
    return 4;
  }

  return 0;
}
