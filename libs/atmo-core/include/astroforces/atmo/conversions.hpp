/**
 * @file conversions.hpp
 * @brief Time and coordinate conversion helpers.
 * @author Watosn
 */
#pragma once

#include <cmath>
#include <ctime>
#include <utility>

#include "astroforces/atmo/constants.hpp"
#include "astroforces/atmo/types.hpp"

namespace astroforces::core {

inline GeodeticPoint spherical_geodetic_from_ecef(const Vec3& ecef_m) {
  constexpr double kPi = 3.1415926535897932384626433832795;
  const double r = norm(ecef_m);
  if (r <= 0.0) {
    return GeodeticPoint{};
  }
  const double lat = std::asin(ecef_m.z / r) * 180.0 / kPi;
  const double lon = std::atan2(ecef_m.y, ecef_m.x) * 180.0 / kPi;
  return GeodeticPoint{.lat_deg = lat, .lon_deg = lon, .alt_m = r - constants::kEarthRadiusWgs84M};
}

inline std::pair<int, double> utc_seconds_to_iyd_sec(double utc_seconds) {
  const std::time_t tt = static_cast<std::time_t>(utc_seconds);
  std::tm tm_utc{};
#if defined(_WIN32)
  gmtime_s(&tm_utc, &tt);
#else
  gmtime_r(&tt, &tm_utc);
#endif
  const int year2 = (tm_utc.tm_year + 1900) % 100;
  const int doy = tm_utc.tm_yday + 1;
  const int iyd = year2 * 1000 + doy;
  const double sec = static_cast<double>(tm_utc.tm_hour * 3600 + tm_utc.tm_min * 60 + tm_utc.tm_sec);
  return {iyd, sec};
}

inline double local_solar_time_hours(double utc_seconds, double lon_deg) {
  const auto iyd_sec = utc_seconds_to_iyd_sec(utc_seconds);
  const double ut_h = iyd_sec.second / 3600.0;
  double stl = std::fmod(ut_h + lon_deg / 15.0, 24.0);
  if (stl < 0.0) {
    stl += 24.0;
  }
  return stl;
}

inline double utc_seconds_to_julian_date_utc(double utc_seconds) {
  return utc_seconds / 86400.0 + 2440587.5;
}

inline double gmst_rad_from_jd_utc(double jd_utc) {
  constexpr double kPi = 3.1415926535897932384626433832795;
  const double t = (jd_utc - 2451545.0) / 36525.0;
  double gmst_deg = 280.46061837 + 360.98564736629 * (jd_utc - 2451545.0) + 0.000387933 * t * t - (t * t * t) / 38710000.0;
  gmst_deg = std::fmod(gmst_deg, 360.0);
  if (gmst_deg < 0.0) {
    gmst_deg += 360.0;
  }
  return gmst_deg * kPi / 180.0;
}

inline Vec3 rotate_z(double angle_rad, const Vec3& v) {
  const double c = std::cos(angle_rad);
  const double s = std::sin(angle_rad);
  return Vec3{
      c * v.x - s * v.y,
      s * v.x + c * v.y,
      v.z,
  };
}

inline Vec3 eci_to_ecef_position(const Vec3& r_eci_m, double utc_seconds) {
  const double gmst = gmst_rad_from_jd_utc(utc_seconds_to_julian_date_utc(utc_seconds));
  return rotate_z(gmst, r_eci_m);
}

inline Vec3 ecef_to_eci_position(const Vec3& r_ecef_m, double utc_seconds) {
  const double gmst = gmst_rad_from_jd_utc(utc_seconds_to_julian_date_utc(utc_seconds));
  return rotate_z(-gmst, r_ecef_m);
}

inline Vec3 eci_to_ecef_velocity(const Vec3& r_eci_m, const Vec3& v_eci_mps, double utc_seconds) {
  const double gmst = gmst_rad_from_jd_utc(utc_seconds_to_julian_date_utc(utc_seconds));
  const Vec3 omega_cross_r_eci{
      -constants::kEarthRotationRateRadPerSec * r_eci_m.y,
      constants::kEarthRotationRateRadPerSec * r_eci_m.x,
      0.0,
  };
  return rotate_z(gmst, v_eci_mps - omega_cross_r_eci);
}

inline Vec3 ecef_to_eci_velocity(const Vec3& r_ecef_m, const Vec3& v_ecef_mps, double utc_seconds) {
  const double gmst = gmst_rad_from_jd_utc(utc_seconds_to_julian_date_utc(utc_seconds));
  const Vec3 r_eci_m = rotate_z(-gmst, r_ecef_m);
  const Vec3 omega_cross_r_eci{
      -constants::kEarthRotationRateRadPerSec * r_eci_m.y,
      constants::kEarthRotationRateRadPerSec * r_eci_m.x,
      0.0,
  };
  return rotate_z(-gmst, v_ecef_mps) + omega_cross_r_eci;
}

}  // namespace astroforces::core
