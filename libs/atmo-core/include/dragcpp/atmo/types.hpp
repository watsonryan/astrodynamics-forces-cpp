/**
 * @file types.hpp
 * @brief Core domain types for drag-cpp.
 * @author Watosn
 */
#pragma once

#include <cmath>
#include <cstdint>

namespace dragcpp::atmo {

enum class Frame : std::uint8_t { ECI, ECEF, NED, BODY };

enum class Status : std::uint8_t { Ok, InvalidInput, NotImplemented, DataUnavailable, NumericalError };

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
};

struct GeodeticPoint {
  double lat_deg{};
  double lon_deg{};
  double alt_m{};
};

struct WeatherIndices {
  double f107{};
  double f107a{};
  double ap{};
  double kp{};
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

}  // namespace dragcpp::atmo
