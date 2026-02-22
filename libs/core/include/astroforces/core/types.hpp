/**
 * @file types.hpp
 * @brief Core domain types for astroforces.
 * @author Watosn
 */
#pragma once

#include <array>
#include <cmath>
#include <cstdint>

namespace astroforces::core {

enum class Frame : std::uint8_t { ECI, ECEF, NED, BODY };

enum class Status : std::uint8_t { Ok, InvalidInput, NotImplemented, DataUnavailable, NumericalError };

enum class WeatherSource : std::uint8_t { Unknown, StaticProvider, CelesTrakLast5YearsCsv };

struct Vec3 {
  double x{};
  double y{};
  double z{};
};

inline Vec3 operator+(const Vec3& a, const Vec3& b) { return Vec3{a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return Vec3{a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator*(double s, const Vec3& v) { return Vec3{s * v.x, s * v.y, s * v.z}; }
inline Vec3 operator*(const Vec3& v, double s) { return s * v; }
inline Vec3 operator/(const Vec3& v, double s) { return Vec3{v.x / s, v.y / s, v.z / s}; }

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double norm(const Vec3& v) { return std::sqrt(dot(v, v)); }

struct Epoch {
  double utc_seconds{};
};

struct StateVector {
  Epoch epoch{};
  Vec3 position_m{};
  Vec3 velocity_mps{};
  Frame frame{Frame::ECI};
  std::array<double, 9> body_from_frame_dcm{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
};

struct GeodeticPoint {
  double lat_deg{};
  double lon_deg{};
  double alt_m{};
};

struct EarthOrientation {
  double xp_rad{};
  double yp_rad{};
  double dut1_s{};
  double lod_s{};
  double dX_rad{};
  double dY_rad{};
};

struct EarthOrientationRate {
  double xp_rad_s{};
  double yp_rad_s{};
  double dut1_s_s{};
  double lod_s_s{};
  double dX_rad_s{};
  double dY_rad_s{};
};

struct CelestialIntermediatePole {
  double x_rad{};
  double y_rad{};
  double s_rad{};
};

struct CelestialIntermediatePoleRate {
  double x_rad_s{};
  double y_rad_s{};
  double s_rad_s{};
};

struct WeatherIndices {
  double f107{};
  double f107a{};
  double ap{};
  double kp{};
  double ap_3h_current{};
  double kp_3h_current{};
  std::array<double, 8> ap_3h_utc{};
  std::array<double, 8> kp_3h_utc{};
  std::array<double, 7> ap_msis_history{};
  bool has_ap_msis_history{};
  bool f107_observed{true};
  bool geomagnetic_observed{true};
  WeatherSource source{WeatherSource::Unknown};
  bool interpolated{};
  bool extrapolated{};
  Status status{Status::Ok};
};

struct AtmosphereSample {
  double density_kg_m3{};
  double temperature_k{};
  Status status{Status::Ok};
};

struct WindSample {
  Vec3 velocity_mps{};
  Frame frame{Frame::ECI};
  Status status{Status::Ok};
};

}  // namespace astroforces::core
