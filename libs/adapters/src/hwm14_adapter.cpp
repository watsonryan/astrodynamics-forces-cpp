/**
 * @file hwm14_adapter.cpp
 * @brief HWM14 adapter implementation.
 * @author Watosn
 */

#include "dragcpp/adapters/hwm14_adapter.hpp"

#include <cmath>
#include <utility>

#include "dragcpp/atmo/conversions.hpp"
#include "hwm14/hwm14.hpp"

namespace astroforces::adapters {

class Hwm14WindAdapter::Impl {
 public:
  explicit Impl(hwm14::Model model) : model_(std::move(model)) {}

  hwm14::Model model_;
};

std::unique_ptr<Hwm14WindAdapter> Hwm14WindAdapter::Create(const Config& config) {
  auto ptr = std::unique_ptr<Hwm14WindAdapter>(new Hwm14WindAdapter(config));
  hwm14::Result<hwm14::Model, hwm14::Error> loaded =
      config.data_dir.empty() ? hwm14::Model::LoadWithSearchPaths() : hwm14::Model::LoadFromDirectory(config.data_dir);
  if (!loaded.has_value()) {
    return ptr;
  }
  ptr->impl_ = std::make_shared<Impl>(loaded.value());
  return ptr;
}

astroforces::atmo::WindSample Hwm14WindAdapter::evaluate(const astroforces::atmo::StateVector& state,
                                                      const astroforces::atmo::WeatherIndices& weather) const {
  if (!impl_) {
    return astroforces::atmo::WindSample{.status = astroforces::atmo::Status::DataUnavailable};
  }
  if (state.frame != astroforces::atmo::Frame::ECEF) {
    return astroforces::atmo::WindSample{.status = astroforces::atmo::Status::InvalidInput};
  }

  const auto geo = astroforces::atmo::spherical_geodetic_from_ecef(state.position_m);
  const auto iyd_sec = astroforces::atmo::utc_seconds_to_iyd_sec(state.epoch.utc_seconds);

  hwm14::Inputs in{};
  in.yyddd = iyd_sec.first;
  in.ut_seconds = iyd_sec.second;
  in.altitude_km = geo.alt_m * 1e-3;
  in.geodetic_lat_deg = geo.lat_deg;
  in.geodetic_lon_deg = geo.lon_deg;
  in.ap3 = weather.ap;

  const auto out = impl_->model_.TotalWinds(in);
  if (!out.has_value()) {
    return astroforces::atmo::WindSample{.status = astroforces::atmo::Status::NumericalError};
  }

  constexpr double kPi = 3.1415926535897932384626433832795;
  const double lat = geo.lat_deg * kPi / 180.0;
  const double lon = geo.lon_deg * kPi / 180.0;
  const double v_n = out.value().meridional_mps;
  const double v_e = out.value().zonal_mps;

  const double vx = -std::sin(lat) * std::cos(lon) * v_n - std::sin(lon) * v_e;
  const double vy = -std::sin(lat) * std::sin(lon) * v_n + std::cos(lon) * v_e;
  const double vz = std::cos(lat) * v_n;

  return astroforces::atmo::WindSample{
      .velocity_mps = astroforces::atmo::Vec3{vx, vy, vz},
      .frame = astroforces::atmo::Frame::ECEF,
      .status = astroforces::atmo::Status::Ok};
}

}  // namespace astroforces::adapters
