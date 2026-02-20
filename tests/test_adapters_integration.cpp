/**
 * @file test_adapters_integration.cpp
 * @brief Adapter-backed integration tests across NRLMSIS, DTM2020, and HWM14.
 * @author Watosn
 */

#include <cmath>
#include <algorithm>
#include <filesystem>
#include <iostream>

#include "dragcpp/adapters/dtm2020_adapter.hpp"
#include "dragcpp/adapters/hwm14_adapter.hpp"
#include "dragcpp/adapters/nrlmsis21_adapter.hpp"
#include "dragcpp/atmo/conversions.hpp"
#include "dragcpp/drag/drag_model.hpp"
#include "dragcpp/sc/spacecraft.hpp"
#include "dragcpp/weather/static_provider.hpp"
#include "msis21/msis21.hpp"

namespace {

bool finite(const dragcpp::atmo::Vec3& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

bool approx_rel(double a, double b, double rel = 1e-12) {
  const double d = std::abs(a - b);
  const double n = std::max(std::abs(b), 1e-30);
  return d / n <= rel;
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

  atmo::WeatherIndices wx_hist = wx;
  wx_hist.ap = 18.0;
  wx_hist.ap_3h_current = 40.0;
  wx_hist.kp_3h_current = 5.0;
  wx_hist.ap_msis_history = {40.0, 36.0, 32.0, 28.0, 22.0, 16.0, 10.0};
  wx_hist.has_ap_msis_history = true;
  const auto nrl_hist_out = nrl->evaluate(state, wx_hist);
  if (nrl_hist_out.status != atmo::Status::Ok) {
    std::cerr << "nrlmsis adapter history evaluation failed\n";
    return 11;
  }

  msis21::Options msis_options{};
  auto msis_model = msis21::Model::load_from_file(nrl_parm, msis_options);
  const auto geo = atmo::spherical_geodetic_from_ecef(state.position_m);
  const auto iyd_sec = atmo::utc_seconds_to_iyd_sec(state.epoch.utc_seconds);
  msis21::Input msis_in{};
  msis_in.iyd = iyd_sec.first;
  msis_in.sec = iyd_sec.second;
  msis_in.alt_km = geo.alt_m * 1.0e-3;
  msis_in.glat_deg = geo.lat_deg;
  msis_in.glon_deg = geo.lon_deg;
  msis_in.stl_hr = atmo::local_solar_time_hours(state.epoch.utc_seconds, geo.lon_deg);
  msis_in.f107a = wx_hist.f107a;
  msis_in.f107 = wx_hist.f107;
  msis_in.ap = wx_hist.ap_3h_current;
  msis_in.ap_history = wx_hist.ap_msis_history;
  msis_in.has_ap_history = true;
  const auto msis_direct_hist = msis_model.evaluate(msis_in);
  if (msis_direct_hist.status != msis21::Status::Ok) {
    std::cerr << "direct nrlmsis history evaluate failed\n";
    return 12;
  }
  if (!approx_rel(nrl_hist_out.density_kg_m3, msis_direct_hist.out.rho * 1000.0) ||
      !approx_rel(nrl_hist_out.temperature_k, msis_direct_hist.out.t)) {
    std::cerr << "nrlmsis adapter parity mismatch for history mode\n";
    return 13;
  }

  msis_in.has_ap_history = false;
  const auto msis_direct_scalar = msis_model.evaluate(msis_in);
  if (msis_direct_scalar.status != msis21::Status::Ok) {
    std::cerr << "direct nrlmsis scalar evaluate failed\n";
    return 14;
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
