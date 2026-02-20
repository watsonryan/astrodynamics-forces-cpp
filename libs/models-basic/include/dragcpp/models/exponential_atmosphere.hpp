/**
 * @file exponential_atmosphere.hpp
 * @brief Basic exponential atmosphere model for integration testing.
 * @author Watosn
 */
#pragma once

#include "dragcpp/atmo/interfaces.hpp"

namespace dragcpp::models {

class ExponentialAtmosphereModel final : public dragcpp::atmo::IAtmosphereModel {
 public:
  ExponentialAtmosphereModel(double rho0_kg_m3, double h0_m, double scale_height_m, double temperature_k)
      : rho0_(rho0_kg_m3), h0_(h0_m), hs_(scale_height_m), t_(temperature_k) {}

  [[nodiscard]] dragcpp::atmo::AtmosphereSample evaluate(const dragcpp::atmo::StateVector& state,
                                                          const dragcpp::atmo::WeatherIndices& weather) const override;

 private:
  double rho0_{};
  double h0_{};
  double hs_{};
  double t_{};
};

class ZeroWindModel final : public dragcpp::atmo::IWindModel {
 public:
  [[nodiscard]] dragcpp::atmo::WindSample evaluate(const dragcpp::atmo::StateVector& state,
                                                    const dragcpp::atmo::WeatherIndices& weather) const override;
};

}  // namespace dragcpp::models
