/**
 * @file nrlmsis21_adapter.hpp
 * @brief Adapter between dragcpp atmosphere interface and NRLMSIS 2.1 model.
 * @author Watosn
 */
#pragma once

#include <filesystem>
#include <memory>

#include "dragcpp/atmo/interfaces.hpp"

namespace msis21 {
class Model;
}

namespace dragcpp::adapters {

class Nrlmsis21AtmosphereAdapter final : public dragcpp::atmo::IAtmosphereModel {
 public:
  struct Config {
    std::filesystem::path parm_file{};
  };

  static std::unique_ptr<Nrlmsis21AtmosphereAdapter> Create(const Config& config);

  [[nodiscard]] dragcpp::atmo::AtmosphereSample evaluate(const dragcpp::atmo::StateVector& state,
                                                          const dragcpp::atmo::WeatherIndices& weather) const override;

 private:
  class Impl;
  explicit Nrlmsis21AtmosphereAdapter(Config config) : config_(std::move(config)) {}

  Config config_{};
  std::shared_ptr<Impl> impl_{};
};

}  // namespace dragcpp::adapters
