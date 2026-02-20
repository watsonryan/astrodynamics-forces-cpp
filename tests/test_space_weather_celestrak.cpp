/**
 * @file test_space_weather_celestrak.cpp
 * @brief CelesTrak CSV weather provider tests.
 * @author Watosn
 */

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "dragcpp/atmo/types.hpp"
#include "dragcpp/weather/celestrak_csv_provider.hpp"

namespace {

bool approx(double a, double b, double tol) { return std::abs(a - b) <= tol; }

bool write_csv(const std::filesystem::path& p) {
  std::ofstream out(p);
  if (!out) {
    return false;
  }
  out << "DATE,BSRN,ND,KP1,KP2,KP3,KP4,KP5,KP6,KP7,KP8,KP_SUM,AP1,AP2,AP3,AP4,AP5,AP6,AP7,AP8,AP_AVG,CP,C9,ISN,F10.7_OBS,F10.7_ADJ,F10.7_DATA_TYPE,F10.7_OBS_CENTER81,F10.7_OBS_LAST81,F10.7_ADJ_CENTER81,F10.7_ADJ_LAST81\n";
  out << "2026-01-01,0,0,10,10,10,10,10,10,10,10,80,6,6,6,6,6,6,6,6,6,0.0,0,0,120.0,120.0,OBS,130.0,130.0,130.0,130.0\n";
  out << "2026-01-02,0,0,20,20,20,20,20,20,20,20,160,14,14,14,14,14,14,14,14,14,0.0,0,0,160.0,160.0,OBS,170.0,170.0,170.0,170.0\n";
  return true;
}

}  // namespace

int main() {
  namespace fs = std::filesystem;
  const auto csv = fs::temp_directory_path() / "dragcpp_celestrak_test.csv";
  if (!write_csv(csv)) {
    std::cerr << "failed to write csv\n";
    return 10;
  }

  const auto provider = dragcpp::weather::CelesTrakCsvSpaceWeatherProvider::Create({.csv_file = csv});

  const dragcpp::atmo::WeatherIndices w0 = provider->at(dragcpp::atmo::Epoch{.utc_seconds = 1767225600.0});  // 2026-01-01
  if (w0.status != dragcpp::atmo::Status::Ok || !approx(w0.f107, 120.0, 1e-12) || !approx(w0.f107a, 130.0, 1e-12) ||
      !approx(w0.ap, 6.0, 1e-12) || !approx(w0.kp, 10.0, 1e-12) || w0.source != dragcpp::atmo::WeatherSource::CelesTrakLast5YearsCsv) {
    std::cerr << "exact day lookup failed\n";
    return 1;
  }

  const dragcpp::atmo::WeatherIndices wm =
      provider->at(dragcpp::atmo::Epoch{.utc_seconds = 1767268800.0});  // 2026-01-01 12:00:00
  if (wm.status != dragcpp::atmo::Status::Ok || !wm.interpolated || wm.extrapolated || !approx(wm.f107, 140.0, 1e-12) ||
      !approx(wm.f107a, 150.0, 1e-12) || !approx(wm.ap, 10.0, 1e-12) || !approx(wm.kp, 15.0, 1e-12)) {
    std::cerr << "interpolation failed\n";
    return 2;
  }

  const dragcpp::atmo::WeatherIndices wb = provider->at(dragcpp::atmo::Epoch{.utc_seconds = 1767139200.0});  // 2025-12-31
  if (wb.status != dragcpp::atmo::Status::Ok || !wb.extrapolated || !approx(wb.f107, 120.0, 1e-12)) {
    std::cerr << "low-side clamp failed\n";
    return 3;
  }

  const dragcpp::atmo::WeatherIndices wa = provider->at(dragcpp::atmo::Epoch{.utc_seconds = 1767312000.0});  // 2026-01-02
  if (wa.status != dragcpp::atmo::Status::Ok || !wa.extrapolated || !approx(wa.f107, 160.0, 1e-12)) {
    std::cerr << "high-side clamp failed\n";
    return 4;
  }

  return 0;
}
