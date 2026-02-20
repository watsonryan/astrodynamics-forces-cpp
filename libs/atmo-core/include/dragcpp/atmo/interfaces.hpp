/**
 * @file interfaces.hpp
 * @brief Core model interfaces for atmosphere/wind/weather.
 * @author Watosn
 */
#pragma once

#include "dragcpp/atmo/types.hpp"

namespace astroforces::atmo {

class ISpaceWeatherProvider {
 public:
  virtual ~ISpaceWeatherProvider() = default;
  [[nodiscard]] virtual WeatherIndices at(const Epoch& epoch) const = 0;
};

class IAtmosphereModel {
 public:
  virtual ~IAtmosphereModel() = default;
  [[nodiscard]] virtual AtmosphereSample evaluate(const StateVector& state,
                                                  const WeatherIndices& weather) const = 0;
};

class IWindModel {
 public:
  virtual ~IWindModel() = default;
  [[nodiscard]] virtual WindSample evaluate(const StateVector& state,
                                            const WeatherIndices& weather) const = 0;
};

}  // namespace astroforces::atmo
