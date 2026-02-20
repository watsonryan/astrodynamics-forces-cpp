/**
 * @file exponential_atmosphere.cpp
 * @brief Basic atmosphere model implementation.
 * @author Watosn
 */

#include "dragcpp/models/exponential_atmosphere.hpp"

#include <cmath>

namespace dragcpp::models {

dragcpp::atmo::AtmosphereSample ExponentialAtmosphereModel::evaluate(const dragcpp::atmo::StateVector& state,
                                                                      const dragcpp::atmo::WeatherIndices& /*weather*/) const {
  const double r = dragcpp::atmo::norm(state.position_m);
  constexpr double kEarthRadiusM = 6378137.0;
  const double alt = r - kEarthRadiusM;
  const double rho = rho0_ * std::exp(-(alt - h0_) / hs_);
  return dragcpp::atmo::AtmosphereSample{.density_kg_m3 = rho, .temperature_k = t_, .status = dragcpp::atmo::Status::Ok};
}

}  // namespace dragcpp::models
