# MAMC: Magnetic Activity Mechanisms Calculator

[![ReadTheDocs](https://img.shields.io/badge/docs-read__the__docs-blue.svg)](https://mamc.readthedocs.io) 
[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![Astropy](https://img.shields.io/badge/powered%20by-Astropy-orange.svg)](https://www.astropy.org/)

**MAMC** is an open-source, object-oriented Python framework designed to evaluate the physical, energetic, and dynamic feasibility of proposed orbital period modulation mechanisms in close eclipsing binary systems—specifically Post-Common Envelope Binaries (PCEBs). 

The software unifies several foundational and modern frameworks to determine if cyclic Eclipse Timing Variations (ETVs) observed in $O-C$ diagrams can be explained by internal magnetic variations, spin-orbit envelope coupling, or non-axisymmetric dynamo mass deformations, or if a third-body (circumbinary planet) hypothesis must be invoked.

---

## 🌌 Implemented Physics Models

MAMC isolates its diagnostics into three independent testing suites:

### 1. The Applegate Mechanism (`ApplegateMechanism`)
Evaluates the energetic feasibility of variations in the active star's magnetic quadrupole moment ($\Delta Q$) based on the numerical benchmarks established by <a href="https://ui.adsabs.harvard.edu/abs/1992ApJ...385..621A/abstract">**Applegate (1992)**</a> based on <a href= "https://ui.adsabs.harvard.edu/abs/2016A%26A...587A..34V/abstract">**Völschow et al. (2016)**]</a> models. It calculates the required energy fraction ($\Delta E / E_{\text{sec}}$) across three separate shell geometries:
* **Thin-Shell Approximation:** Following [Tian et al. (2009)](https://ui.adsabs.harvard.edu/abs/2009Ap%26SS.319..119T/abstract), testing energy constraints based on absolute period variations.
* **Constant-Density Model:** A global structural approximation assuming a uniform density profile.
* **Two-Zone Model:** A realistic treatment modeling an inner dense core and an outer fluid convective envelope exchanging angular momentum. It tracks the dimensionless structural Applegate parameter ($A$), where physical solutions strictly require $A \le 1.0$.

### 2. Lanza Spin-Orbit Envelope Coupling (`LanzaMechanism`)
Evaluates the radiant energy constraints and mechanical shearing indicators for the convective shell velocity modulations proposed by [**Lanza (2020)**](https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.1820L/abstract). The framework maps observed timing amplitudes down to cyclical core-shell rotation fluctuations ($\Delta \Omega / \Omega$) and checks the total mechanical work against the active companion's total radiant luminosity budget ($L_{\text{sec}}$).

### 3. Azimuthal Dynamo Wave Model (`AzimuthalDynamoWave`)
Implements the analytical scaling framework described by [**Navarrete et al. (2026)**](https://ui.adsabs.harvard.edu/abs/2026arXiv260427609N/abstract) to test if non-axisymmetric quadrupole moment changes—driven by strong, azimuthally migrating magnetic fields in rapidly rotating convective stars ($\alpha^2$ dynamo regime)—can natively reproduce observed ETV amplitudes. It maps principal inertia tensor variations ($\Delta I_{xx}, \Delta I_{yy}, \Delta I_{zz}$) under rapid rotation to test prolate geometry configurations relative to the rotational axis.

---

## 📂 Repository Architecture

The project is structured modularly to completely decouple the computational physics equations from target data sets and pipeline execution scripts:

```text
├── mamc.py            # Computational Core: Contains physics classes and conversion logic.
├── target_params.py   # Target Catalog: Structured @dataclass containing system parameters.
└── example_mamc.py    # Pipeline Runner: Main execution script generating separated diagnostic reports.
