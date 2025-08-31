

<details>
<summary><h1>two_plate_mgrf: Modern Python Package for Electrical Double Layer Simulations (Two Plates)</h1></summary>

<p align="center">
  <img src="https://img.shields.io/badge/python-3.8%2B-blue" alt="Python Version">
  <img src="https://img.shields.io/badge/license-MIT-green" alt="License">
  <img src="https://img.shields.io/badge/build-passing-brightgreen" alt="Build Status">
</p>

---

## Overview

**two_plate_mgrf** is a Python package for simulating the electrical double layer (EDL) structure between two uniformly charged plates using the modified Gaussian renormalized fluctuation (MGRF) theory. The code is built on top of the [Dedalus](https://github.com/DedalusProject/dedalus) spectral PDE solver and implements advanced iterative and parallel algorithms for high-accuracy EDL modeling.

This package enables researchers to:
- Solve the MGRF equations for planar EDLs with arbitrary salt mixtures in two-plate geometry
- Reproduce published results from recent theoretical works
- Extend and adapt the code for new physical scenarios
- This version is for systems with uniform dielectric permittivity, check "main" branch for  code for systems with dielectric contrast.

<details>
<summary>References & Citations</summary>

If you use this code, please cite:

- Agrawal & Wang, [Phys. Rev. Lett. 2022, 129, 228001](https://doi.org/10.1103/PhysRevLett.129.228001)
- Agrawal & Wang, [J. Chem. Theory Comput. 2022, 18, 6271–6280](https://doi.org/10.1021/acs.jctc.2c00607)
- Nikhil R. Agrawal, Ravtej Kaur, Carlo Carraro, and Rui Wang [arXiv:2306.10137](https://doi.org/10.48550/arXiv.2306.10137)

</details>

---

## Features

- Full solution of the MGRF theory for planar EDLs (two-plate geometry)
- Support for arbitrary salt mixtures and dielectric contrasts
- Parallelized Green’s function and self-energy calculations
- Modular, extensible codebase with clear separation of physics and numerics
- Output in HDF5 format for easy post-processing
- Reproducibility: all parameters and results are saved for each run

---

## Installation

1. **Clone the repository:**
	```bash
	git clone https://github.com/nikhil0165/two_plates_mgrf.git
	cd two_plates_mgrf
	```
2. **Set up a Python environment (recommended):**
	```bash
	conda create -n mgrf python=3.8
	conda activate mgrf
	pip install -r requirements.txt
	# Or install dependencies manually (see packages.py)
	```
3. **Install Dedalus:**
	Follow [Dedalus installation instructions](https://dedalus-project.readthedocs.io/en/latest/).

---

## Usage

You can run the main simulation scripts as follows:

```bash
# Run with a new initial guess (PB solution)
python simulator_pb.py physical_param.py

# Run using a saved MGRF solution as initial guess
python simulator.py physical_param.py
```

All output is saved in the `results*` folders as HDF5 files, with all parameters and profiles included.

**Note:**
`simulator_pb.py` saves the solution of the modified Gaussian renormalized fluctuation theory in a .h5 file for the input parameters given in the two *_param.py files, using the PB solution as the initial guess. `simulator.py` uses a saved solution of `mgrf_2plate.py` as the initial guess. The physical variables for this saved solution and the final parameters for which you want the double-layer structure can be set separately using the file `physical_param.py`. The variables deciding which saved solution to choose as initial guess end with `_in_d`, for example: `sigma_in_d`.

---

## File Structure & Module Guide

<details>
<summary>Click to expand file/module descriptions</summary>

- **numerical_param.py**: Numerical parameters (grid size, tolerances, mixing ratios, etc.) for all solvers.
- **physical_param.py**: Physical system parameters (concentrations, valencies, radii, dielectric constants, etc.).
- **dh_2plate.py**: Linearized Debye–Hückel solver for initial guess.
- **pb_2plate.py**: Nonlinear Poisson–Boltzmann solver for mean-field EDL structure.
- **mgrf_2plate.py**: Main MGRF solver for 2-plate geometry.
- **num_concn.py**: Functions for computing ion concentration profiles and Jacobians.
- **selfe_2plate.py**: Self-energy calculations for the interface (2-plate).
- **selfe_bulk.py**: Self-energy calculations for the bulk.
- **greens_function_2plate.py**: Green’s function (Fourier transform) for the interface.
- **greens_function_bulk.py**: Green’s function for the bulk.
- **calculate.py**: Utility functions for screening length, ionic strength, charge density, interpolation, etc.
- **energy_2plate.py**: Grand free energy calculations for interface and bulk.
- **simulator_pb.py**: Script to solve and save MGRF solution using PB as initial guess.
- **simulator.py**: Script to solve and save MGRF solution using a previous MGRF result as initial guess.
- **packages.py**: Centralized imports for all required Python libraries.

</details>

---

## Contributing

Contributions, bug reports, and feature requests are welcome! Please open an issue or submit a pull request.

---

## Contact

Developed by Nikhil Agrawal in the lab of Prof. Rui Wang, Pitzer Center for Theoretical Chemistry, University of California, Berkeley, USA.

For questions, suggestions, or collaboration, contact: <nikhilagrawal0165@gmail.com>

</details>

