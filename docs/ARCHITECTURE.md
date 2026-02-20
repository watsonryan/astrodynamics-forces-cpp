# astrodynamics-forces-cpp Architecture

## Modules
- `atmo-core`: canonical types and interfaces
- `space-weather`: providers for F10.7/Ap/Kp and related inputs
- `models-basic`: temporary baseline models for integration tests
- `adapters`: integration seam for NRLMSIS/DTM/HWM wrappers
- `sc-props`: spacecraft macro geometry and drag-relevant properties
- `drag-core`: drag acceleration pipeline + generic perturbation interfaces
- `apps/drag-cli`: executable entrypoint

## Force Abstraction
- `dragcpp::forces::IPerturbationModel`: common force contribution interface.
- `dragcpp::forces::PerturbationStack`: additive combiner for all perturbation models.
- `dragcpp::drag::DragPerturbationModel`: drag implementation of the generic interface.

## Design Rules
- Preserve model kernels when integrating external repos.
- Make units and frames explicit in all interfaces.
- Keep adapters thin; keep model-specific assumptions local.
