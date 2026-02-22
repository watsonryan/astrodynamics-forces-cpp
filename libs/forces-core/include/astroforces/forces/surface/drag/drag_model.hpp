/**
 * @file drag_model.hpp
 * @brief Drag acceleration pipeline.
 * @author Watosn
 */
#pragma once

#include "astroforces/core/interfaces.hpp"
#include "astroforces/core/eop.hpp"
#include "astroforces/core/cip.hpp"
#include "astroforces/sc/spacecraft.hpp"

namespace astroforces::forces {

enum class DragFrameTransformMode : unsigned char {
  ApproxGmst,
  StrictGcrfItrf,
};

struct DragResult {
  astroforces::core::Vec3 acceleration_mps2{};
  astroforces::core::Vec3 relative_velocity_mps{};
  double density_kg_m3{};
  double temperature_k{};
  double relative_speed_mps{};
  double dynamic_pressure_pa{};
  double area_m2{};
  double cd{};
  astroforces::core::WeatherIndices weather{};
  astroforces::core::Status status{astroforces::core::Status::Ok};
};

class DragAccelerationModel {
 public:
  DragAccelerationModel(const astroforces::core::ISpaceWeatherProvider& weather,
                        const astroforces::core::IAtmosphereModel& atmosphere,
                        const astroforces::core::IWindModel& wind,
                        DragFrameTransformMode transform_mode = DragFrameTransformMode::ApproxGmst,
                        const astroforces::core::eop::Series* eop_series = nullptr,
                        const astroforces::core::cip::Series* cip_series = nullptr)
      : weather_(weather),
        atmosphere_(atmosphere),
        wind_(wind),
        transform_mode_(transform_mode),
        eop_series_(eop_series),
        cip_series_(cip_series) {}

  [[nodiscard]] DragResult evaluate(const astroforces::core::StateVector& state,
                                    const astroforces::sc::SpacecraftProperties& sc) const;

 private:
  const astroforces::core::ISpaceWeatherProvider& weather_;
  const astroforces::core::IAtmosphereModel& atmosphere_;
  const astroforces::core::IWindModel& wind_;
  DragFrameTransformMode transform_mode_{};
  const astroforces::core::eop::Series* eop_series_{nullptr};
  const astroforces::core::cip::Series* cip_series_{nullptr};
};

}  // namespace astroforces::forces
