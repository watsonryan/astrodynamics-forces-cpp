/**
 * @file drag_model.cpp
 * @brief Drag acceleration pipeline implementation.
 * @author Watosn
 */

#include "astroforces/forces/surface/drag/drag_model.hpp"

#include <cmath>

#include <Eigen/Dense>

#include "astroforces/core/transforms.hpp"
#include "astroforces/forces/surface/surface_force.hpp"

namespace astroforces::forces {

DragResult DragAccelerationModel::evaluate(const astroforces::core::StateVector& state,
                                           const astroforces::sc::SpacecraftProperties& sc) const {
  if (sc.mass_kg <= 0.0) {
    return DragResult{.status = astroforces::core::Status::InvalidInput};
  }
  if (state.frame != astroforces::core::Frame::ECI && state.frame != astroforces::core::Frame::ECEF) {
    return DragResult{.status = astroforces::core::Status::InvalidInput};
  }

  const auto w = weather_.at(state.epoch);
  if (w.status != astroforces::core::Status::Ok) {
    return DragResult{.status = w.status};
  }

  astroforces::core::StateVector eval_state = state;
  if (state.frame == astroforces::core::Frame::ECI) {
    eval_state.frame = astroforces::core::Frame::ECEF;
    eval_state.position_m = astroforces::core::eci_to_ecef_position(state.position_m, state.epoch.utc_seconds);
    eval_state.velocity_mps =
        astroforces::core::eci_to_ecef_velocity(state.position_m, state.velocity_mps, state.epoch.utc_seconds);
  }

  const auto a = atmosphere_.evaluate(eval_state, w);
  if (a.status != astroforces::core::Status::Ok || a.density_kg_m3 < 0.0) {
    return DragResult{.status = a.status};
  }

  const auto wind = wind_.evaluate(eval_state, w);
  if (wind.status != astroforces::core::Status::Ok) {
    return DragResult{.status = wind.status};
  }

  astroforces::core::Vec3 wind_state_mps{};
  if (state.frame == astroforces::core::Frame::ECEF) {
    if (wind.frame == astroforces::core::Frame::ECEF) {
      wind_state_mps = wind.velocity_mps;
    } else if (wind.frame == astroforces::core::Frame::ECI) {
      const auto r_eci = astroforces::core::ecef_to_eci_position(eval_state.position_m, state.epoch.utc_seconds);
      wind_state_mps = astroforces::core::eci_to_ecef_velocity(r_eci, wind.velocity_mps, state.epoch.utc_seconds);
    } else {
      return DragResult{.status = astroforces::core::Status::InvalidInput};
    }
  } else {
    if (wind.frame == astroforces::core::Frame::ECI) {
      wind_state_mps = wind.velocity_mps;
    } else if (wind.frame == astroforces::core::Frame::ECEF) {
      wind_state_mps =
          astroforces::core::ecef_to_eci_velocity(eval_state.position_m, wind.velocity_mps, state.epoch.utc_seconds);
    } else {
      return DragResult{.status = astroforces::core::Status::InvalidInput};
    }
  }

  const auto vrel = state.velocity_mps - wind_state_mps;
  double speed = 0.0;
  const auto flow_dir_frame = astroforces::forces::unit_direction(vrel, &speed);
  if (!std::isfinite(speed)) {
    return DragResult{.status = astroforces::core::Status::NumericalError};
  }

  astroforces::core::Vec3 flow_dir_body{};
  if (speed > 0.0) {
    const Eigen::Vector3d flow_frame(flow_dir_frame.x, flow_dir_frame.y, flow_dir_frame.z);
    const Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> body_from_frame(state.body_from_frame_dcm.data());
    const Eigen::Vector3d flow_body = body_from_frame * flow_frame;
    flow_dir_body = astroforces::core::Vec3{flow_body.x(), flow_body.y(), flow_body.z()};
  }

  const double q_pa = 0.5 * a.density_kg_m3 * speed * speed;
  const auto sf = astroforces::forces::evaluate_surface_force(sc,
                                                              flow_dir_frame,
                                                              flow_dir_body,
                                                              q_pa,
                                                              sc.cd,
                                                              astroforces::sc::SurfaceCoeffModel::Drag,
                                                              -1.0);
  if (sf.status != astroforces::core::Status::Ok) {
    return DragResult{.status = sf.status};
  }

  return DragResult{.acceleration_mps2 = sf.acceleration_mps2,
                    .relative_velocity_mps = vrel,
                    .density_kg_m3 = a.density_kg_m3,
                    .temperature_k = a.temperature_k,
                    .relative_speed_mps = speed,
                    .dynamic_pressure_pa = q_pa,
                    .area_m2 = sf.area_m2,
                    .cd = sf.coeff,
                    .weather = w,
                    .status = astroforces::core::Status::Ok};
}

}  // namespace astroforces::forces
