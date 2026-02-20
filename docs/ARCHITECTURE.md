# Architecture

## Modules
- `atmo-core`: canonical types and interfaces
- `space-weather`: providers for F10.7/Ap/Kp and related inputs
- `models-basic`: temporary baseline models for integration tests
- `adapters`: integration seam for NRLMSIS/DTM/HWM wrappers
- `sc-props`: spacecraft macro geometry and drag-relevant properties
- `drag-core`: drag acceleration physics pipeline
- `apps/drag-cli`: executable entrypoint

## Design Rules
- Preserve model kernels when integrating external repos.
- Make units and frames explicit in all interfaces.
- Keep adapters thin; keep model-specific assumptions local.
