/**
 * @file exponential_atmosphere.cpp
 * @brief Basic atmosphere model implementation.
 * @author Watosn
 */

#include "astroforces/models/exponential_atmosphere.hpp"

#include <cmath>

#include "astroforces/atmo/constants.hpp"

namespace astroforces::models {

astroforces::atmo::AtmosphereSample ExponentialAtmosphereModel::evaluate(const astroforces::atmo::StateVector& state,
                                                                      const astroforces::atmo::WeatherIndices& /*weather*/) const {
  const double r = astroforces::atmo::norm(state.position_m);
  const double alt = r - astroforces::atmo::constants::kEarthRadiusWgs84M;
  const double rho = rho0_ * std::exp(-(alt - h0_) / hs_);
  return astroforces::atmo::AtmosphereSample{.density_kg_m3 = rho, .temperature_k = t_, .status = astroforces::atmo::Status::Ok};
}

}  // namespace astroforces::models
