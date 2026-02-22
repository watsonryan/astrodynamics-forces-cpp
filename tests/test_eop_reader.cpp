/**
 * @file test_eop_reader.cpp
 * @brief EOP finals parser/interpolation tests.
 * @author Watosn
 */

#include <cmath>
#include <filesystem>

#include <spdlog/spdlog.h>

#include "astroforces/core/eop.hpp"

namespace {

bool approx(double a, double b, double tol) { return std::abs(a - b) <= tol; }

}  // namespace

int main() {
  using namespace astroforces::core;
  using namespace astroforces::core::eop;

#ifndef ASTROFORCES_SOURCE_DIR
  spdlog::error("ASTROFORCES_SOURCE_DIR missing");
  return 10;
#endif

  const std::filesystem::path sample = std::filesystem::path(ASTROFORCES_SOURCE_DIR) / "tests" / "data" / "eop_finals_sample.txt";
  const auto series = Series::load_iers_finals(sample);
  if (series.empty() || series.size() != 3) {
    spdlog::error("unexpected EOP series size");
    return 1;
  }

  const auto first = series.at_mjd_utc(60676.0);
  if (!first.has_value()) {
    spdlog::error("missing first EOP row");
    return 2;
  }
  if (!approx(first->xp_rad / astroforces::core::constants::kArcsecToRad, 0.123456, 1e-12)
      || !approx(first->yp_rad / astroforces::core::constants::kArcsecToRad, 0.234567, 1e-12)
      || !approx(first->dut1_s, 0.1, 1e-12)
      || !approx(first->lod_s, 0.0009, 1e-15)) {
    spdlog::error("first row values mismatch");
    return 3;
  }

  const auto mid = series.at_mjd_utc(60676.5);
  if (!mid.has_value()) {
    spdlog::error("missing interpolated row");
    return 4;
  }
  if (!approx(mid->xp_rad / astroforces::core::constants::kArcsecToRad, 0.173456, 1e-9)
      || !approx(mid->yp_rad / astroforces::core::constants::kArcsecToRad, 0.284567, 1e-9)
      || !approx(mid->dut1_s, 0.11, 1e-12)
      || !approx(mid->lod_s, 0.0010, 1e-15)) {
    spdlog::error("interpolation mismatch");
    return 5;
  }

  const auto clamped = series.at_mjd_utc(60690.0);
  if (!clamped.has_value() || !approx(clamped->xp_rad / astroforces::core::constants::kArcsecToRad, 0.323456, 1e-12)) {
    spdlog::error("clamp-to-end mismatch");
    return 6;
  }

  return 0;
}
