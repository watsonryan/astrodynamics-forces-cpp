# drag-cpp Master Plan

## 1. Goal
Build a unified C++ drag modeling platform that computes spacecraft drag acceleration from:
- Orbital state and epoch
- Spacecraft geometry + surface/material properties
- Atmospheric density model(s): NRLMSIS 2.1, DTM2020
- Wind model(s): HWM14 (and later optional models)
- Space weather inputs (F10.7, F10.7a, Ap/Kp, etc.)

Success criteria:
- Reproducible, validated drag acceleration outputs
- Preserved numerical parity for imported models
- Clean, extensible architecture for adding models and datasets

## 2. Architecture Strategy
Use a monorepo with strict library boundaries. Do not rewrite validated kernels unless required.

### Proposed layout
- `libs/atmo-core`
  - Shared domain types: epoch, frames, geodetic/geocentric, units
  - Common interfaces for atmosphere/wind/weather models
- `libs/space-weather`
  - Readers/parsers/cache/interpolation for indices and forecast data
- `libs/nrlmsis21`
  - Imported/embedded from `nrlmsis-2_1` with minimal API adapter layer
- `libs/dtm2020`
  - Imported/embedded from `dtm2020` with adapter layer
- `libs/hwm14`
  - Imported/embedded from `hwm14` with adapter layer
- `libs/drag-core`
  - Relative velocity, Cd model hooks, area projection, acceleration equation
- `libs/sc-props`
  - Spacecraft geometry, macro-model surfaces, material coefficients
- `apps/drag-cli`
  - Batch/single-point runs, profile generation, diagnostics
- `tests/`
  - Unit, model-parity, integration, regression, performance
- `docs/`
  - Validation reports, architecture, data assumptions, API guides

## 3. Canonical Interfaces (define first)
Define stable interfaces before migration.

### `AtmosphereModel`
Input:
- Epoch
- Position (ECEF/ECI + frame tag)
- Optional precomputed geodetic fields
- Space weather snapshot

Output:
- Total density (kg/m^3)
- Species (optional)
- Temperature(s)
- Status/diagnostics

### `WindModel`
Output:
- Neutral wind vector in a declared frame (m/s)

### `SpaceWeatherProvider`
Input:
- Epoch
Output:
- Struct with all required indices + quality flags + provenance

### `DragModel`
Input:
- State vector
- Spacecraft properties
- Atmosphere + wind outputs
Output:
- Drag acceleration vector (m/s^2)
- Component breakdown for debug

## 4. Migration Principles
1. Preserve model kernels first; wrap second.
2. No hidden unit conversions: all conversions explicit and tested.
3. Every adapter has parity tests against source repo baselines.
4. Keep model-specific assumptions local to each adapter.
5. Prefer deterministic execution paths and reproducible builds.

## 5. Phased Plan

### Phase 0: Bootstrap
- Initialize repo structure + top-level CMake presets aligned with your standard.
- Add shared lint/format/test tooling.
- Add baseline CI skeleton (build + unit tests only initially).

Deliverables:
- Buildable empty modules
- Standard presets and toolchain docs

### Phase 1: Core Contracts and Units
- Implement shared types:
  - `Epoch`, `StateVector`, `FrameTag`, `GeodeticCoord`, `WeatherIndices`
- Implement unit-safe wrappers or strict naming conventions.
- Add frame/unit conversion utilities with golden tests.

Deliverables:
- `atmo-core` API stable
- Unit/frame test suite

### Phase 2: Import Existing Repos as Libraries
- Bring in `nrlmsis-2_1`, `dtm2020`, `hwm14` as module libraries.
- Add thin adapters implementing the common interfaces.
- Keep original tests passing inside each module.

Deliverables:
- Model adapters compile and run
- Existing model repos’ tests preserved

### Phase 3: Space Weather Layer
- Implement `space-weather` provider:
  - File-based readers (historical, operational)
  - Interpolation/windowing rules
  - Caching strategy
- Add provenance and quality flags for each value.

Deliverables:
- Unified `WeatherIndices` source for all models
- Reader/interpolation regression tests

### Phase 4: Drag Physics Core
- Implement drag equation pipeline:
  - Relative velocity = spacecraft inertial velocity transformed into atmosphere frame minus winds
  - Dynamic pressure terms
  - Area projection from geometry/macro-surface model
  - Cd model plugin architecture
- Start with fixed Cd mode, then optional advanced Cd models.

Deliverables:
- `drag-core` acceleration API
- Verified unit tests on analytic toy cases

### Phase 5: End-to-End Integration
- Build `drag-cli`:
  - Single-point evaluation
  - Time-series batch run
  - CSV/JSON output
- Create scenarios to compare NRLMSIS vs DTM + optional HWM contribution.

Deliverables:
- Operational end-to-end executable
- Reproducible scenario scripts

### Phase 6: Validation and Parity Lock
- Model parity gates:
  - Keep existing per-model parity checks
- Cross-model sanity checks:
  - Physical ranges
  - Continuity across altitude/time windows
- Drag regression suite:
  - Golden trajectories or stepwise accelerations

Deliverables:
- Validation report in `docs/validation.md`
- CI gates for parity + integration

### Phase 7: Performance and Production Hardening
- Profile hotspots in adapter + drag pipeline.
- Add benchmark targets and performance budgets.
- Add robust error taxonomy and structured logging.

Deliverables:
- Perf baselines and thresholds
- Stable release candidate

## 6. Testing Matrix

### Unit tests
- Units, frames, math primitives, interpolation, geometry projection

### Parity tests
- NRLMSIS adapter parity against existing `nrlmsis-2_1` baseline
- DTM parity against `dtm2020` baseline
- HWM parity against `hwm14` baseline

### Integration tests
- Atmosphere + wind + weather + drag acceleration end-to-end

### Regression tests
- Golden cases for mission-like states and weather windows

### Performance tests
- Throughput and latency benchmarks for single-point and batch modes

## 7. Key Risks and Mitigations

### Risk: unit/frame mismatch
Mitigation:
- Strongly typed APIs or explicit suffix naming
- Mandatory conversion tests and frame tags

### Risk: model assumption mismatch
Mitigation:
- Adapter-local normalization logic
- Model-specific docs and precondition checks

### Risk: silent numerical drift
Mitigation:
- Preserve operation order in migrated kernels
- Tight parity tests and per-compiler baseline monitoring

### Risk: dependency/data sprawl
Mitigation:
- Centralized data loading layer + immutable weather snapshots

## 8. Suggested Milestones
- M1 (1-2 weeks): Phases 0-1 complete
- M2 (2-4 weeks): Phases 2-3 complete
- M3 (2-3 weeks): Phase 4 complete
- M4 (2 weeks): Phase 5 complete
- M5 (ongoing): Phases 6-7 hardening

## 9. Definition of Done
- All imported model parity tests pass at agreed tolerances.
- End-to-end drag acceleration outputs validated on representative scenarios.
- Performance baselines documented and gated.
- Clear API and operations docs for users.

## 10. Immediate Next Actions
1. Create repo skeleton and top-level CMake targets.
2. Copy in model repos as isolated libs (or submodules) with adapters.
3. Freeze interface contracts in `atmo-core`.
4. Add first end-to-end drag acceleration path using one density model + fixed Cd.
