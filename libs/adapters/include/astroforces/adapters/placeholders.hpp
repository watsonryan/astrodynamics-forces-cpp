/**
 * @file placeholders.hpp
 * @brief Placeholder adapters for model integration.
 * @author Watosn
 */
#pragma once

namespace astroforces::adapters {

// Integration points for imported external model repos:
// - nrlmsis-2_1
// - dtm2020
// - hwm14
// Real adapters should implement astroforces::atmo interfaces.

inline constexpr const char* kAdapterStatus = "placeholders_ready";

}  // namespace astroforces::adapters
