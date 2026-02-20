# drag-cpp

Unified C++20 drag modeling platform for spacecraft drag acceleration.

## Scope
- Shared atmospheric/wind/weather interfaces
- Drag acceleration core (relative velocity + ballistic term)
- Spacecraft geometry/surface property hooks
- Adapter integration points for:
  - NRLMSIS 2.1 (`nrlmsis-2_1`)
  - DTM2020 (`dtm2020`)
  - HWM14 (`hwm14`)

## Build
```bash
cmake --preset macos-debug
cmake --build --preset macos-debug
ctest --preset macos-debug --output-on-failure
```

## CLI
```bash
./build/macos-debug/drag_cli 6778137 0 0 0 7670 0 1000000000
```

## Next Integration Steps
1. Replace `models-basic` with adapter-backed model bundle wiring.
2. Implement real space weather readers/interpolation in `libs/space-weather`.
3. Add frame/attitude handling and surface-resolved drag in `libs/drag-core`.
