/**
 * @file celestrak_csv_provider.cpp
 * @brief CelesTrak SW-Last5Years CSV provider implementation.
 * @author Watosn
 */

#include "dragcpp/weather/celestrak_csv_provider.hpp"

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>

namespace dragcpp::weather {
namespace {

constexpr std::size_t kMinCelesTrakColumns = 31;
constexpr std::size_t kDateCol = 0;
constexpr std::size_t kKpSumCol = 11;
constexpr std::size_t kApAvgCol = 20;
constexpr std::size_t kF107ObsCol = 24;
constexpr std::size_t kF107ObsCenter81Col = 27;

int days_from_civil(int y, unsigned m, unsigned d) {
  y -= static_cast<int>(m <= 2);
  const int era = (y >= 0 ? y : y - 399) / 400;
  const unsigned yoe = static_cast<unsigned>(y - era * 400);
  const unsigned doy = (153U * (m + (m > 2 ? -3U : 9U)) + 2U) / 5U + d - 1U;
  const unsigned doe = yoe * 365U + yoe / 4U - yoe / 100U + doy;
  return era * 146097 + static_cast<int>(doe) - 719468;
}

double ymd_to_utc_seconds(int y, unsigned m, unsigned d) {
  return static_cast<double>(days_from_civil(y, m, d)) * 86400.0;
}

std::vector<std::string> split_csv_line(const std::string& line) {
  std::vector<std::string> fields;
  fields.reserve(40);
  std::string token;
  std::stringstream ss(line);
  while (std::getline(ss, token, ',')) {
    fields.push_back(token);
  }
  return fields;
}

bool parse_date(const std::string& text, int& year, unsigned& month, unsigned& day) {
  if (text.size() != 10 || text[4] != '-' || text[7] != '-') {
    return false;
  }
  try {
    year = std::stoi(text.substr(0, 4));
    month = static_cast<unsigned>(std::stoul(text.substr(5, 2)));
    day = static_cast<unsigned>(std::stoul(text.substr(8, 2)));
  } catch (...) {
    return false;
  }
  return month >= 1U && month <= 12U && day >= 1U && day <= 31U;
}

bool parse_double(const std::string& text, double& value) {
  try {
    value = std::stod(text);
  } catch (...) {
    return false;
  }
  return true;
}

dragcpp::atmo::WeatherIndices make_indices(const CelesTrakCsvSpaceWeatherProvider::DailySample& s, bool interpolated,
                                           bool extrapolated) {
  return dragcpp::atmo::WeatherIndices{
      .f107 = s.f107_obs,
      .f107a = s.f107_obs_center81,
      .ap = s.ap_avg,
      .kp = s.kp_avg,
      .source = dragcpp::atmo::WeatherSource::CelesTrakLast5YearsCsv,
      .interpolated = interpolated,
      .extrapolated = extrapolated,
      .status = dragcpp::atmo::Status::Ok};
}

}  // namespace

std::unique_ptr<CelesTrakCsvSpaceWeatherProvider> CelesTrakCsvSpaceWeatherProvider::Create(const Config& config) {
  std::ifstream in(config.csv_file);
  if (!in) {
    return std::unique_ptr<CelesTrakCsvSpaceWeatherProvider>(
        new CelesTrakCsvSpaceWeatherProvider(std::vector<DailySample>{}));
  }

  std::vector<DailySample> samples;
  std::string line;
  bool header_consumed = false;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    if (!header_consumed) {
      header_consumed = true;
      if (line.find("DATE") != std::string::npos) {
        continue;
      }
    }

    const auto fields = split_csv_line(line);
    if (fields.size() < kMinCelesTrakColumns) {
      continue;
    }

    int year = 0;
    unsigned month = 0;
    unsigned day = 0;
    if (!parse_date(fields[kDateCol], year, month, day)) {
      continue;
    }

    double kp_sum = 0.0;
    double ap_avg = 0.0;
    double f107 = 0.0;
    double f107a = 0.0;
    if (!parse_double(fields[kKpSumCol], kp_sum) || !parse_double(fields[kApAvgCol], ap_avg) ||
        !parse_double(fields[kF107ObsCol], f107) || !parse_double(fields[kF107ObsCenter81Col], f107a)) {
      continue;
    }

    samples.push_back(DailySample{
        .day_start_utc_s = ymd_to_utc_seconds(year, month, day),
        .f107_obs = f107,
        .f107_obs_center81 = f107a,
        .ap_avg = ap_avg,
        .kp_avg = kp_sum / 8.0,
    });
  }

  std::sort(samples.begin(), samples.end(),
            [](const DailySample& a, const DailySample& b) { return a.day_start_utc_s < b.day_start_utc_s; });
  return std::unique_ptr<CelesTrakCsvSpaceWeatherProvider>(new CelesTrakCsvSpaceWeatherProvider(std::move(samples)));
}

dragcpp::atmo::WeatherIndices CelesTrakCsvSpaceWeatherProvider::at(const dragcpp::atmo::Epoch& epoch) const {
  if (samples_.empty()) {
    return dragcpp::atmo::WeatherIndices{
        .source = dragcpp::atmo::WeatherSource::CelesTrakLast5YearsCsv,
        .status = dragcpp::atmo::Status::DataUnavailable};
  }

  const double t = epoch.utc_seconds;
  if (t <= samples_.front().day_start_utc_s) {
    return make_indices(samples_.front(), false, true);
  }
  if (t >= samples_.back().day_start_utc_s) {
    return make_indices(samples_.back(), false, true);
  }

  const auto it = std::upper_bound(samples_.begin(), samples_.end(), t,
                                   [](double ts, const DailySample& s) { return ts < s.day_start_utc_s; });
  if (it == samples_.begin()) {
    return make_indices(samples_.front(), false, true);
  }

  const auto& hi = *it;
  const auto& lo = *(it - 1);
  const double dt = hi.day_start_utc_s - lo.day_start_utc_s;
  if (dt <= 0.0) {
    return make_indices(lo, false, false);
  }
  const double alpha = std::clamp((t - lo.day_start_utc_s) / dt, 0.0, 1.0);

  DailySample interp{};
  interp.f107_obs = lo.f107_obs + alpha * (hi.f107_obs - lo.f107_obs);
  interp.f107_obs_center81 = lo.f107_obs_center81 + alpha * (hi.f107_obs_center81 - lo.f107_obs_center81);
  interp.ap_avg = lo.ap_avg + alpha * (hi.ap_avg - lo.ap_avg);
  interp.kp_avg = lo.kp_avg + alpha * (hi.kp_avg - lo.kp_avg);

  return make_indices(interp, alpha > 0.0 && alpha < 1.0, false);
}

}  // namespace dragcpp::weather
