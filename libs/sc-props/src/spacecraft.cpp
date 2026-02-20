/**
 * @file spacecraft.cpp
 * @brief Spacecraft geometry helper implementation.
 * @author Watosn
 */

#include "dragcpp/sc/spacecraft.hpp"

#include <algorithm>

namespace dragcpp::sc {

double projected_area_m2(const SpacecraftProperties& sc, const dragcpp::atmo::Vec3& flow_dir_body) {
  if (!sc.use_surface_model || sc.surfaces.empty()) {
    return sc.reference_area_m2;
  }

  double area = 0.0;
  for (const auto& s : sc.surfaces) {
    const double c = -dragcpp::atmo::dot(s.normal_body, flow_dir_body);
    area += s.area_m2 * std::max(0.0, c);
  }
  return area;
}

}  // namespace dragcpp::sc
