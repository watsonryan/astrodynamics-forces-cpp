/**
 * @file hwm14_adapter.hpp
 * @brief Adapter between dragcpp wind interface and HWM14 model.
 * @author Watosn
 */
#pragma once

#include <filesystem>
#include <memory>

#include "dragcpp/atmo/interfaces.hpp"

namespace hwm14 {
class Model;
}

namespace dragcpp::adapters {

class Hwm14WindAdapter final : public dragcpp::atmo::IWindModel {
 public:
  struct Config {
    std::filesystem::path data_dir{};
  };

  static std::unique_ptr<Hwm14WindAdapter> Create(const Config& config);

  [[nodiscard]] dragcpp::atmo::WindSample evaluate(const dragcpp::atmo::StateVector& state,
                                                    const dragcpp::atmo::WeatherIndices& weather) const override;

 private:
  class Impl;
  explicit Hwm14WindAdapter(Config config) : config_(std::move(config)) {}

  Config config_{};
  std::shared_ptr<Impl> impl_{};
};

}  // namespace dragcpp::adapters
