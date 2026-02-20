/**
 * @file dtm2020_adapter.hpp
 * @brief Adapter between dragcpp atmosphere interface and DTM2020 model.
 * @author Watosn
 */
#pragma once

#include <filesystem>
#include <memory>

#include "dragcpp/atmo/interfaces.hpp"

namespace dtm2020 {
class Dtm2020Operational;
}

namespace dragcpp::adapters {

class Dtm2020AtmosphereAdapter final : public dragcpp::atmo::IAtmosphereModel {
 public:
  struct Config {
    std::filesystem::path coeff_file{};
  };

  static std::unique_ptr<Dtm2020AtmosphereAdapter> Create(const Config& config);

  [[nodiscard]] dragcpp::atmo::AtmosphereSample evaluate(const dragcpp::atmo::StateVector& state,
                                                          const dragcpp::atmo::WeatherIndices& weather) const override;

 private:
  class Impl;
  explicit Dtm2020AtmosphereAdapter(Config config) : config_(std::move(config)) {}

  Config config_{};
  std::shared_ptr<Impl> impl_{};
};

}  // namespace dragcpp::adapters
