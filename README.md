# mamc
# Magnetic Activity Mechanism Calculator

This repository provides a python-based analytical framework implemented in Jupyter Notebook format to explore the physical and energetic feasibility of magnetic mechanisms in explaining Eclipse Timing Variations (ETVs) observed in Post-Common-Envelope Binaries (PCEBs).

---

## 🚀 Overview

The codebase evaluates whether observed orbital period modulations in close binary systems can be driven by a secondary star's magnetic cycles. It provides implementations for two core theoretical frameworks:

1. **The Applegate Mechanism**: Calculates the energy ratio ($\Delta E / E_{\mathrm{sec}}$) required to drive orbital variations across three distinct structural approximations (Völschow et al. 2016).
2. **The Azimuthal Dynamo Wave (ADW) Model**: Estimates theoretical $O-C$ timing semi-amplitudes resulting from non-axisymmetric quadrupole moment variations (Navarrete et al. 2026).

---

## 📋 Core Features

* **Multiple Applegate Formulations**: Includes implementations for the Thin-Shell approximation (Tian et al. 2009), the Constant-Density model, and the Two-Zone model (Völschow et al. 2016)[cite: 1].
* **ADW Scaling Analysis**: Maps inertia tensor variations from 3D MHD simulations to non-axisymmetric quadrupole moments to approximate geometric timing shifts under tidal synchronization assumptions[cite: 1].
* **Literature Benchmarking**: Contains pre-loaded parameters and verification checks for several well-studied systems[cite: 1]:
  * **NN Ser**, **HW Vir**, and **QS Vir** (reproducing Völschow et al. 2016)[cite: 1]
  * **NY Vir** (reproducing Esmer et al. 2023)[cite: 1]
  * **DD CrB** (reproducing Baştürk et al. 2026)[cite: 1]
* **Robust Conversions**: Built natively around `astropy.constants` and `astropy.units` to handle seamless, error-free CGS conversions[cite: 1].

---

## 🛠️ Getting Started

### Prerequisites
Make sure you have a Python environment with the following dependencies installed:
```bash
pip install numpy astropy notebook
