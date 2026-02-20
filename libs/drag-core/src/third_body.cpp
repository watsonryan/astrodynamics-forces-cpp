/**
 * @file third_body.cpp
 * @brief Third-body perturbation model implementation.
 * @author Watosn
 */

#include "dragcpp/forces/third_body.hpp"

#include <cmath>

#include "dragcpp/atmo/conversions.hpp"
#include "jpl_eph/jpl_eph.hpp"

namespace astroforces::forces {
namespace {

astroforces::atmo::Vec3 to_vec3(const std::array<double, 6>& pv) {
  return astroforces::atmo::Vec3{pv[0], pv[1], pv[2]};
}

astroforces::atmo::Status map_jpl_error(const jpl::eph::Status& s) {
  switch (s.code) {
    case jpl::eph::ErrorCode::kInvalidArgument:
      return astroforces::atmo::Status::InvalidInput;
    case jpl::eph::ErrorCode::kIo:
    case jpl::eph::ErrorCode::kCorruptFile:
    case jpl::eph::ErrorCode::kOutOfRange:
    case jpl::eph::ErrorCode::kUnsupported:
      return astroforces::atmo::Status::DataUnavailable;
    case jpl::eph::ErrorCode::kOk:
    default:
      return astroforces::atmo::Status::NumericalError;
  }
}

astroforces::atmo::Vec3 third_body_direct_indirect(const astroforces::atmo::Vec3& r_sc_m,
                                               const astroforces::atmo::Vec3& r_tb_m,
                                               double mu_m3_s2) {
  const astroforces::atmo::Vec3 rho = r_tb_m - r_sc_m;
  const double rho_norm = astroforces::atmo::norm(rho);
  const double rtb_norm = astroforces::atmo::norm(r_tb_m);
  if (rho_norm <= 0.0 || rtb_norm <= 0.0) {
    return astroforces::atmo::Vec3{};
  }
  const double rho3 = rho_norm * rho_norm * rho_norm;
  const double rtb3 = rtb_norm * rtb_norm * rtb_norm;
  return (mu_m3_s2 / rho3) * rho - (mu_m3_s2 / rtb3) * r_tb_m;
}

}  // namespace

ThirdBodyPerturbationModel::ThirdBodyPerturbationModel(Config config) : config_(std::move(config)) {}

std::unique_ptr<ThirdBodyPerturbationModel> ThirdBodyPerturbationModel::Create(const Config& config) {
  auto out = std::unique_ptr<ThirdBodyPerturbationModel>(new ThirdBodyPerturbationModel(config));
  auto opened = jpl::eph::Ephemeris::Open(config.ephemeris_file.string());
  if (!opened.has_value()) {
    return out;
  }
  out->ephemeris_ = opened.value();
  out->workspace_ = std::make_shared<jpl::eph::Workspace>();
  return out;
}

PerturbationContribution ThirdBodyPerturbationModel::evaluate(const PerturbationRequest& request) const {
  PerturbationContribution out{};
  out.name = config_.name;
  out.type = PerturbationType::ThirdBody;

  if (!ephemeris_ || !workspace_) {
    out.status = astroforces::atmo::Status::DataUnavailable;
    return out;
  }
  if (request.state.frame != astroforces::atmo::Frame::ECI) {
    out.status = astroforces::atmo::Status::InvalidInput;
    return out;
  }

  // Approximate UTC as TDB for initial force-model integration.
  const double jed_tdb = astroforces::atmo::utc_seconds_to_julian_date_utc(request.state.epoch.utc_seconds);

  if (config_.use_sun) {
    const auto sun = ephemeris_->PlephSi(jed_tdb, jpl::eph::Body::Sun, jpl::eph::Body::Earth, false, *workspace_);
    if (!sun.has_value()) {
      out.status = map_jpl_error(sun.error());
      return out;
    }
    out.acceleration_mps2 =
        out.acceleration_mps2 + third_body_direct_indirect(request.state.position_m, to_vec3(sun.value().pv), config_.mu_sun_m3_s2);
  }

  if (config_.use_moon) {
    const auto moon = ephemeris_->PlephSi(jed_tdb, jpl::eph::Body::Moon, jpl::eph::Body::Earth, false, *workspace_);
    if (!moon.has_value()) {
      out.status = map_jpl_error(moon.error());
      return out;
    }
    out.acceleration_mps2 =
        out.acceleration_mps2 + third_body_direct_indirect(request.state.position_m, to_vec3(moon.value().pv), config_.mu_moon_m3_s2);
  }

  out.status = astroforces::atmo::Status::Ok;
  return out;
}

}  // namespace astroforces::forces
