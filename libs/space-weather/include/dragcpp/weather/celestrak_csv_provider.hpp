/**
 * @file celestrak_csv_provider.hpp
 * @brief Space weather provider backed by CelesTrak SW-Last5Years CSV.
 * @author Watosn
 */
#pragma once

#include <filesystem>
#include <memory>
#include <vector>

#include "dragcpp/atmo/interfaces.hpp"

namespace dragcpp::weather {

class CelesTrakCsvSpaceWeatherProvider final : public dragcpp::atmo::ISpaceWeatherProvider {
 public:
  struct DailySample {
    double day_start_utc_s{};
    double f107_obs{};
    double f107_obs_center81{};
    double ap_avg{};
    double kp_avg{};
  };

  struct Config {
    std::filesystem::path csv_file{};
  };

  static std::unique_ptr<CelesTrakCsvSpaceWeatherProvider> Create(const Config& config);

  [[nodiscard]] dragcpp::atmo::WeatherIndices at(const dragcpp::atmo::Epoch& epoch) const override;

 private:
  explicit CelesTrakCsvSpaceWeatherProvider(std::vector<DailySample> samples) : samples_(std::move(samples)) {}

  std::vector<DailySample> samples_{};
};

}  // namespace dragcpp::weather
