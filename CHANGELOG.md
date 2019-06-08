# Change Log

## [0.1.2] - 2019.06.08

### Changed

- Added an American-spelt version of the BravaisLattice attribute `centring_type`: `centering_type`.
- Lattice parameters `a`, `b`, `c`, `alpha`/`α`, `beta`/`β` and `gamma`/`γ`, `lattice_system` and `centring_type`/`centering_type` are now properties and therefore cannot be modified after instantiation.

### Fixed

- `row_or_column` attribute is now respected by `lattice_sites` and `lattice_sites_frac` attributes.
- Deprecation warning from `pyyaml` has been resolved.

## [0.1.1] - 2019.06.03

### Changed

- This is the initial "complete" version; PyPI "development status" has been changed from "3 - Alpha" to "4 - Beta".
