/**
 * @file spacecraft.cpp
 * @brief Spacecraft geometry helper implementation.
 * @author Watosn
 */

#include "dragcpp/sc/spacecraft.hpp"

#include <algorithm>

namespace astroforces::sc {

AeroProjection projected_area_and_cd(const SpacecraftProperties& sc, const astroforces::atmo::Vec3& flow_dir_body) {
  if (!sc.use_surface_model || sc.surfaces.empty()) {
    return AeroProjection{.area_m2 = sc.reference_area_m2, .cd_effective = sc.cd};
  }

  double area = 0.0;
  double weighted_cd = 0.0;
  for (const auto& s : sc.surfaces) {
    const double c = -astroforces::atmo::dot(s.normal_body, flow_dir_body);
    const double proj = s.area_m2 * std::max(0.0, c);
    area += proj;
    const double cd_base = (s.cd > 0.0) ? s.cd : sc.cd;
    double cd_surface = cd_base;
    if (s.specularity > 0.0 || s.accommodation > 0.0) {
      const double incidence = std::clamp(c, 0.0, 1.0);
      const double spec = std::clamp(s.specularity, 0.0, 1.0);
      const double accom = std::clamp(s.accommodation, 0.0, 1.0);
      const double modifier = 1.0 + 0.25 * accom * (1.0 + incidence) - 0.15 * spec * (1.0 - incidence);
      cd_surface = cd_base * modifier;
    }
    weighted_cd += cd_surface * proj;
  }
  const double cd_effective = (area > 0.0) ? (weighted_cd / area) : sc.cd;
  return AeroProjection{.area_m2 = area, .cd_effective = cd_effective};
}

double projected_area_m2(const SpacecraftProperties& sc, const astroforces::atmo::Vec3& flow_dir_body) {
  return projected_area_and_cd(sc, flow_dir_body).area_m2;
}

}  // namespace astroforces::sc
