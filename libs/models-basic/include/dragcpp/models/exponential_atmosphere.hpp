/**
 * @file exponential_atmosphere.hpp
 * @brief Basic exponential atmosphere model for integration testing.
 * @author Watosn
 */
#pragma once

#include "dragcpp/atmo/interfaces.hpp"

namespace astroforces::models {

class ExponentialAtmosphereModel final : public astroforces::atmo::IAtmosphereModel {
 public:
  ExponentialAtmosphereModel(double rho0_kg_m3, double h0_m, double scale_height_m, double temperature_k)
      : rho0_(rho0_kg_m3), h0_(h0_m), hs_(scale_height_m), t_(temperature_k) {}

  [[nodiscard]] astroforces::atmo::AtmosphereSample evaluate(const astroforces::atmo::StateVector& state,
                                                          const astroforces::atmo::WeatherIndices& weather) const override;

 private:
  double rho0_{};
  double h0_{};
  double hs_{};
  double t_{};
};

class ZeroWindModel final : public astroforces::atmo::IWindModel {
 public:
  [[nodiscard]] astroforces::atmo::WindSample evaluate(const astroforces::atmo::StateVector& state,
                                                    const astroforces::atmo::WeatherIndices& weather) const override;
};

}  // namespace astroforces::models
