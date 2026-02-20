/**
 * @file static_provider.hpp
 * @brief Constant-value weather provider.
 * @author Watosn
 */
#pragma once

#include "dragcpp/atmo/interfaces.hpp"

namespace dragcpp::weather {

class StaticSpaceWeatherProvider final : public dragcpp::atmo::ISpaceWeatherProvider {
 public:
  explicit StaticSpaceWeatherProvider(dragcpp::atmo::WeatherIndices indices) : indices_(indices) {
    indices_.source = dragcpp::atmo::WeatherSource::StaticProvider;
  }
  [[nodiscard]] dragcpp::atmo::WeatherIndices at(const dragcpp::atmo::Epoch& /*epoch*/) const override {
    return indices_;
  }

 private:
  dragcpp::atmo::WeatherIndices indices_{};
};

}  // namespace dragcpp::weather
