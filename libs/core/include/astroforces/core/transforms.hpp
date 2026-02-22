/**
 * @file transforms.hpp
 * @brief Shared time and coordinate transform helpers.
 * @author Watosn
 */
#pragma once

#include <array>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <utility>

#include "astroforces/core/constants.hpp"
#include "astroforces/core/types.hpp"
#include "astroforces/core/math_utils.hpp"
#include "astroforces/core/sofa_utils.hpp"

namespace astroforces::core {

struct RotationWithDerivative {
  Mat3 r{};
  Mat3 dr{};
};

inline Vec3 cross(const Vec3& a, const Vec3& b) {
  return vec_cross(a, b);
}

inline GeodeticPoint spherical_geodetic_from_ecef(const Vec3& ecef_m) {
  const double r = norm(ecef_m);
  if (r <= 0.0) {
    return GeodeticPoint{};
  }
  const double lat = std::asin(ecef_m.z / r) * 180.0 / constants::kPi;
  const double lon = std::atan2(ecef_m.y, ecef_m.x) * 180.0 / constants::kPi;
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
  return utc_seconds / constants::kSecondsPerDay + 2440587.5;
}

inline double tai_minus_utc_seconds(const double utc_seconds) {
  // Leap-second epochs (UTC) and resulting TAI-UTC seconds.
  constexpr std::array<std::pair<double, double>, 28> kLeapTable{{
      {63072000.0, 10.0},     // 1972-01-01
      {78796800.0, 11.0},     // 1972-07-01
      {94694400.0, 12.0},     // 1973-01-01
      {126230400.0, 13.0},    // 1974-01-01
      {157766400.0, 14.0},    // 1975-01-01
      {189302400.0, 15.0},    // 1976-01-01
      {220924800.0, 16.0},    // 1977-01-01
      {252460800.0, 17.0},    // 1978-01-01
      {283996800.0, 18.0},    // 1979-01-01
      {315532800.0, 19.0},    // 1980-01-01
      {362793600.0, 20.0},    // 1981-07-01
      {394329600.0, 21.0},    // 1982-07-01
      {425865600.0, 22.0},    // 1983-07-01
      {489024000.0, 23.0},    // 1985-07-01
      {567993600.0, 24.0},    // 1988-01-01
      {631152000.0, 25.0},    // 1990-01-01
      {662688000.0, 26.0},    // 1991-01-01
      {709948800.0, 27.0},    // 1992-07-01
      {741484800.0, 28.0},    // 1993-07-01
      {773020800.0, 29.0},    // 1994-07-01
      {820454400.0, 30.0},    // 1996-01-01
      {867715200.0, 31.0},    // 1997-07-01
      {915148800.0, 32.0},    // 1999-01-01
      {1136073600.0, 33.0},   // 2006-01-01
      {1230768000.0, 34.0},   // 2009-01-01
      {1341100800.0, 35.0},   // 2012-07-01
      {1435708800.0, 36.0},   // 2015-07-01
      {1483228800.0, 37.0},   // 2017-01-01
  }};

  double tai_utc = 10.0;
  for (const auto& [epoch_utc_s, delta] : kLeapTable) {
    if (utc_seconds >= epoch_utc_s) {
      tai_utc = delta;
    } else {
      break;
    }
  }
  return tai_utc;
}

inline double tt_minus_utc_seconds(const double utc_seconds) {
  return tai_minus_utc_seconds(utc_seconds) + 32.184;
}

inline double utc_seconds_to_julian_date_tt(const double utc_seconds) {
  return utc_seconds_to_julian_date_utc(utc_seconds) + tt_minus_utc_seconds(utc_seconds) / constants::kSecondsPerDay;
}

inline double earth_rotation_angle_rad(double jd_ut1) {
  return sofa::era00(jd_ut1);
}

inline double gmst_rad_from_jd_utc(double jd_utc) {
  const double t = (jd_utc - constants::kJ2000Jd) / 36525.0;
  double gmst_deg = 280.46061837 + 360.98564736629 * (jd_utc - constants::kJ2000Jd) + 0.000387933 * t * t
                    - (t * t * t) / 38710000.0;
  gmst_deg = std::fmod(gmst_deg, 360.0);
  if (gmst_deg < 0.0) {
    gmst_deg += 360.0;
  }
  return gmst_deg * constants::kDegToRad;
}

inline Vec3 rotate_z(double angle_rad, const Vec3& v) {
  return mat_vec(sofa::rot_z(angle_rad), v);
}

// Approximate GMST-only transform (no precession/nutation/polar motion). Use
// gcrf_to_itrf_* with EOP/CIP inputs for production-quality navigation.
inline Vec3 eci_to_ecef_position(const Vec3& r_eci_m, double utc_seconds) {
  const double gmst = gmst_rad_from_jd_utc(utc_seconds_to_julian_date_utc(utc_seconds));
  return rotate_z(gmst, r_eci_m);
}

// Approximate GMST-only transform (no precession/nutation/polar motion). Use
// itrf_to_gcrf_* with EOP/CIP inputs for production-quality navigation.
inline Vec3 ecef_to_eci_position(const Vec3& r_ecef_m, double utc_seconds) {
  const double gmst = gmst_rad_from_jd_utc(utc_seconds_to_julian_date_utc(utc_seconds));
  return rotate_z(-gmst, r_ecef_m);
}

inline Vec3 eci_to_ecef_velocity(const Vec3& r_eci_m, const Vec3& v_eci_mps, double utc_seconds) {
  const double gmst = gmst_rad_from_jd_utc(utc_seconds_to_julian_date_utc(utc_seconds));
  const Vec3 omega_cross_r_eci{-constants::kEarthRotationRateRadPerSec * r_eci_m.y,
                                constants::kEarthRotationRateRadPerSec * r_eci_m.x,
                                0.0};
  return rotate_z(gmst, v_eci_mps - omega_cross_r_eci);
}

inline Vec3 ecef_to_eci_velocity(const Vec3& r_ecef_m, const Vec3& v_ecef_mps, double utc_seconds) {
  const double gmst = gmst_rad_from_jd_utc(utc_seconds_to_julian_date_utc(utc_seconds));
  const Vec3 r_eci_m = rotate_z(-gmst, r_ecef_m);
  const Vec3 omega_cross_r_eci{-constants::kEarthRotationRateRadPerSec * r_eci_m.y,
                                constants::kEarthRotationRateRadPerSec * r_eci_m.x,
                                0.0};
  return rotate_z(-gmst, v_ecef_mps) + omega_cross_r_eci;
}

inline double tio_locator_sp_rad(const double jd_tt) {
  return sofa::sp00(jd_tt);
}

inline Mat3 c2i_from_xys(const CelestialIntermediatePole& cip) {
  return sofa::c2ixys(cip.x_rad, cip.y_rad, cip.s_rad);
}

inline Mat3 polar_motion_matrix(const EarthOrientation& eop, const double sp_rad) {
  return sofa::pom00(eop.xp_rad, eop.yp_rad, sp_rad);
}

inline Mat3 mat_add(const Mat3& a, const Mat3& b) {
  Mat3 c{};
  for (int i = 0; i < 9; ++i) {
    c.v[static_cast<std::size_t>(i)] = a.v[static_cast<std::size_t>(i)] + b.v[static_cast<std::size_t>(i)];
  }
  return c;
}

inline Mat3 mat_scale(double s, const Mat3& a) {
  Mat3 c{};
  for (int i = 0; i < 9; ++i) {
    c.v[static_cast<std::size_t>(i)] = s * a.v[static_cast<std::size_t>(i)];
  }
  return c;
}

inline double tio_locator_sp_rate_rad_s() {
  return (-47e-6 * constants::kArcsecToRad) / (36525.0 * constants::kSecondsPerDay);
}

inline Mat3 gcrf_to_itrf_rotation(
    const double jd_utc,
    const double jd_tt,
    const CelestialIntermediatePole& cip,
    const EarthOrientation& eop) {
  const double jd_ut1 = jd_utc + eop.dut1_s / constants::kSecondsPerDay;
  CelestialIntermediatePole corrected = cip;
  corrected.x_rad += eop.dX_rad;
  corrected.y_rad += eop.dY_rad;
  const Mat3 rc2i = c2i_from_xys(corrected);
  const double era = earth_rotation_angle_rad(jd_ut1);
  const Mat3 rpom = polar_motion_matrix(eop, tio_locator_sp_rad(jd_tt));
  return sofa::c2tcio(rc2i, era, rpom);
}

inline RotationWithDerivative gcrf_to_itrf_rotation_with_derivative_exact(
    const double jd_utc,
    const double jd_tt,
    const CelestialIntermediatePole& cip,
    const CelestialIntermediatePoleRate& cip_rate,
    const EarthOrientation& eop,
    const EarthOrientationRate& eop_rate) {
  const double x = cip.x_rad + eop.dX_rad;
  const double y = cip.y_rad + eop.dY_rad;
  const double x_dot = cip_rate.x_rad_s + eop_rate.dX_rad_s;
  const double y_dot = cip_rate.y_rad_s + eop_rate.dY_rad_s;
  const double s_cio = cip.s_rad;
  const double s_dot = cip_rate.s_rad_s;
  const bool cip_rate_zero = (x_dot == 0.0) && (y_dot == 0.0) && (s_dot == 0.0);
  const bool pm_rate_zero = (eop_rate.xp_rad_s == 0.0) && (eop_rate.yp_rad_s == 0.0);

  const double r2 = x * x + y * y;
  const double r = (r2 > 0.0) ? std::sqrt(r2) : 0.0;
  const double e = (r2 > 0.0) ? std::atan2(y, x) : 0.0;
  const double e_dot = (r2 > 0.0) ? ((x * y_dot - y * x_dot) / r2) : 0.0;
  const double r_dot = (r > 0.0) ? ((x * x_dot + y * y_dot) / r) : 0.0;
  const double d = (r > 0.0) ? std::asin(r) : 0.0;
  const double d_dot = (r > 0.0) ? (r_dot / std::sqrt(std::max(1e-30, 1.0 - r2))) : 0.0;

  const double a1 = -(e + s_cio);
  const double a2 = d;
  const double a3 = e;
  const double a1_dot = -(e_dot + s_dot);
  const double a2_dot = d_dot;
  const double a3_dot = e_dot;

  const Mat3 rz1 = sofa::rot_z(a1);
  const Mat3 ry2 = sofa::rot_y(a2);
  const Mat3 rz3 = sofa::rot_z(a3);
  const Mat3 drz1 = mat_scale(a1_dot, sofa::drot_z_da(a1));
  const Mat3 dry2 = mat_scale(a2_dot, sofa::drot_y_da(a2));
  const Mat3 drz3 = mat_scale(a3_dot, sofa::drot_z_da(a3));

  const Mat3 ry2_rz3 = mat_mul(ry2, rz3);
  const Mat3 rc2i = mat_mul(rz1, ry2_rz3);
  Mat3 drc2i{};
  if (!cip_rate_zero) {
    const Mat3 t_drc2i_1 = mat_mul(drz1, ry2_rz3);
    const Mat3 t_drc2i_2 = mat_mul(rz1, mat_mul(dry2, rz3));
    const Mat3 t_drc2i_3 = mat_mul(rz1, mat_mul(ry2, drz3));
    drc2i = mat_add(t_drc2i_1, mat_add(t_drc2i_2, t_drc2i_3));
  }

  const double jd_ut1 = jd_utc + eop.dut1_s / constants::kSecondsPerDay;
  const double era = earth_rotation_angle_rad(jd_ut1);
  const double omega = constants::kEarthRotationRateRadPerSec / (1.0 + eop.lod_s / constants::kSecondsPerDay);
  const Mat3 r3era = sofa::rot_z(era);
  const Mat3 dr3era = mat_scale(omega, sofa::drot_z_da(era));

  const double sp = tio_locator_sp_rad(jd_tt);
  const double sp_dot = tio_locator_sp_rate_rad_s();
  const Mat3 rzp = sofa::rot_z(sp);
  const Mat3 ryp = sofa::rot_y(-eop.xp_rad);
  const Mat3 rxp = sofa::rot_x(-eop.yp_rad);
  const Mat3 drzp = mat_scale(sp_dot, sofa::drot_z_da(sp));
  const Mat3 dryp = mat_scale(-eop_rate.xp_rad_s, sofa::drot_y_da(-eop.xp_rad));
  const Mat3 drxp = mat_scale(-eop_rate.yp_rad_s, sofa::drot_x_da(-eop.yp_rad));

  const Mat3 ryp_rxp = mat_mul(ryp, rxp);
  const Mat3 rpom = mat_mul(rzp, ryp_rxp);
  const Mat3 t_drpom_1 = mat_mul(drzp, ryp_rxp);
  Mat3 drpom = t_drpom_1;
  if (!pm_rate_zero) {
    const Mat3 t_drpom_2 = mat_mul(rzp, mat_mul(dryp, rxp));
    const Mat3 t_drpom_3 = mat_mul(rzp, mat_mul(ryp, drxp));
    drpom = mat_add(t_drpom_1, mat_add(t_drpom_2, t_drpom_3));
  }

  RotationWithDerivative out{};
  const Mat3 r3era_rc2i = mat_mul(r3era, rc2i);
  out.r = mat_mul(rpom, r3era_rc2i);
  const Mat3 t_dr_1 = mat_mul(drpom, r3era_rc2i);
  const Mat3 t_dr_2 = mat_mul(rpom, mat_mul(dr3era, rc2i));
  out.dr = mat_add(t_dr_1, t_dr_2);
  if (!cip_rate_zero) {
    const Mat3 t_dr_3 = mat_mul(rpom, mat_mul(r3era, drc2i));
    out.dr = mat_add(out.dr, t_dr_3);
  }
  return out;
}

inline RotationWithDerivative gcrf_to_itrf_rotation_with_derivative(
    const double jd_utc,
    const double jd_tt,
    const CelestialIntermediatePole& cip,
    const EarthOrientation& eop) {
  return gcrf_to_itrf_rotation_with_derivative_exact(jd_utc, jd_tt, cip, CelestialIntermediatePoleRate{}, eop, EarthOrientationRate{});
}

inline Vec3 gcrf_to_itrf_position(
    const Vec3& r_gcrf,
    const double jd_utc,
    const double jd_tt,
    const CelestialIntermediatePole& cip,
    const EarthOrientation& eop) {
  return mat_vec(gcrf_to_itrf_rotation_with_derivative(jd_utc, jd_tt, cip, eop).r, r_gcrf);
}

inline Vec3 itrf_to_gcrf_position(
    const Vec3& r_itrf,
    const double jd_utc,
    const double jd_tt,
    const CelestialIntermediatePole& cip,
    const EarthOrientation& eop) {
  const Mat3 rt = mat_transpose(gcrf_to_itrf_rotation_with_derivative(jd_utc, jd_tt, cip, eop).r);
  return mat_vec(rt, r_itrf);
}

inline Vec3 gcrf_to_itrf_velocity(
    const Vec3& r_gcrf,
    const Vec3& v_gcrf,
    const double jd_utc,
    const double jd_tt,
    const CelestialIntermediatePole& cip,
    const EarthOrientation& eop) {
  const RotationWithDerivative rd = gcrf_to_itrf_rotation_with_derivative(jd_utc, jd_tt, cip, eop);
  return mat_vec(rd.r, v_gcrf) + mat_vec(rd.dr, r_gcrf);
}

inline Vec3 itrf_to_gcrf_velocity(
    const Vec3& r_itrf,
    const Vec3& v_itrf,
    const double jd_utc,
    const double jd_tt,
    const CelestialIntermediatePole& cip,
    const EarthOrientation& eop) {
  const RotationWithDerivative rd = gcrf_to_itrf_rotation_with_derivative(jd_utc, jd_tt, cip, eop);
  const Mat3 rt = mat_transpose(rd.r);
  const Vec3 r_gcrf = mat_vec(rt, r_itrf);
  const Vec3 correction = mat_vec(rd.dr, r_gcrf);
  return mat_vec(rt, v_itrf - correction);
}

}  // namespace astroforces::core
