/**
 * @file zero_wind.cpp
 * @brief Zero neutral wind model implementation.
 * @author Watosn
 */

#include "dragcpp/models/exponential_atmosphere.hpp"

namespace astroforces::models {

astroforces::atmo::WindSample ZeroWindModel::evaluate(const astroforces::atmo::StateVector& state,
                                                   const astroforces::atmo::WeatherIndices& /*weather*/) const {
  return astroforces::atmo::WindSample{
      .velocity_mps = astroforces::atmo::Vec3{},
      .frame = state.frame,
      .status = astroforces::atmo::Status::Ok};
}

}  // namespace astroforces::models
