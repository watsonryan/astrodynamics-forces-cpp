/**
 * @file zero_wind.cpp
 * @brief Zero neutral wind model implementation.
 * @author Watosn
 */

#include "dragcpp/models/exponential_atmosphere.hpp"

namespace dragcpp::models {

dragcpp::atmo::WindSample ZeroWindModel::evaluate(const dragcpp::atmo::StateVector& state,
                                                   const dragcpp::atmo::WeatherIndices& /*weather*/) const {
  return dragcpp::atmo::WindSample{
      .velocity_mps = dragcpp::atmo::Vec3{},
      .frame = state.frame,
      .status = dragcpp::atmo::Status::Ok};
}

}  // namespace dragcpp::models
