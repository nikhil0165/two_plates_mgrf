<summary><h1>two_plate_mgrf: Modern Python Package for Electrical Double Layer Simulations in Slit Geometries</h1></summary>

<p align="center">
  <img src="https://img.shields.io/badge/python-3.8%2B-blue" alt="Python Version">
  <img src="https://img.shields.io/badge/license-MIT-green" alt="License">
  <img src="https://img.shields.io/badge/build-passing-brightgreen" alt="Build Status">
</p>

---

## Overview

**two_plate_mgrf** is a Python package for simulating the electrical double layer (EDL) structure in slit geometries formed by two uniformly charged parallel plates, using the modified Gaussian renormalized fluctuation (MGRF) theory. The code is built on top of the [Dedalus](https://github.com/DedalusProject/dedalus) spectral PDE solver and implements advanced iterative and parallel algorithms for high-accuracy EDL modeling.

This package enables researchers to:
- Solve the MGRF equations for confined EDLs in slit pores with arbitrary salt mixtures and dielectric environments
- Study double-layer interactions and confinement effects beyond mean-field theory
- Reproduce published results from recent theoretical works
- Extend and adapt the code for new physical scenarios

<details>
<summary>References & Citations</summary>

The code implements equations derived in:
- Agrawal and Wang, [Phys. Rev. Lett. 2022, 129, 228001](https://doi.org/10.1103/PhysRevLett.129.228001)
- Agrawal and Wang, [J. Chem. Theory Comput. 2022, 18, 6271–6280](https://doi.org/10.1021/acs.jctc.2c00607)

The iteration scheme for solving non-linear equations is partially adapted from:
- Xu and Maggs, [J. Comp. Phys. 275 (2014): 310-322](https://doi.org/10.1016/j.jcp.2014.07.004)

This code can be used to reproduce results presented in
- Nikhil R. Agrawal, Ravtej Kaur, Carlo Carraro, and Rui Wang [arXiv:2306.10137](https://doi.org/10.48550/arXiv.2306.10137)
- Nikhil R. Agrawal, Carlo Carraro, and Rui Wang [J. Chem. Phys. 161, 204902 (2024)](https://doi.org/10.1063/5.0235611)

</details>

---

## Features

- Full solution of the MGRF theory for slit (two-plate) geometries
- Support for arbitrary salt mixtures and dielectric contrasts
- Parallelized Green’s function and self-energy calculations
- Modular, extensible codebase with clear separation of physics and numerics
- Output in HDF5 format for easy post-processing
- Reproducibility: all parameters and results are saved for each run

---

## Installation

1. **Clone the repository:**
	```bash
	git clone https://github.com/nikhil0165/two_plate_mgrf.git
	cd two_plate_mgrf
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




## File Descriptions

### Input Files
- **`numerical_param.py`**: Specifies numerical parameters like grid points, tolerance criteria, mixing ratios, and NCC cutoffs for Dedalus.
- **`physical_param.py`**: Defines physical parameters such as salt concentrations, ion valencies, Born radii, domain size, dielectric constants, and derived non-dimensional variables.

### Core Solvers
- **`dh_2plate.py`**: Solves the linearized mean-field Poisson-Boltzmann (Debye-Hückel) theory for EDLs.
- **`pb_2plate.py`**: Solves the full mean-field Poisson-Boltzmann theory for EDLs.
- **`mgrf_2plate.py`**: Computes the solution to the MGRF theory for a two-plate system.

### Supporting Modules
- **`num_concn.py`**: Functions for calculating concentration profiles and coefficients for iterative schemes.
- **`selfe_2plate.py`**: Calculates self-energy profiles for a two-plate system.
- **`selfe_bulk.py`**: Computes self-energy for bulk solutions.
- **`greens_function_2plate.py` & `greens_function_bulk.py`**: Evaluate Fourier transforms of Green's functions in the interface and bulk.
- **`calculate.py`**: Utility functions for properties like screening length, ionic strength, and charge density profiles.
- **`energy_2plate.py`**: Calculates the grand free energy for both mean-field PB and MGRF theories.

### Simulators
- **`simulator_pb.py`**: Saves the MGRF solution using `pb_2plate.py` as the initial guess.
- **`simulator.py`**: Saves the MGRF solution using a previously saved solution as the initial guess.

### Utilities
- **`packages.py`**: Contains all required Python library imports.

## How to Run
Run the code using one of the following commands:
```bash
python simulator.py physical_param.py
python simulator_pb.py physical_param.py
```

## Dependencies
Ensure you have the required libraries installed. We recommend creating a dedicated Conda environment for this package.

## Contact
Developed by **Nikhil Agrawal** in the lab of Prof. Rui Wang, Pitzer Center for Theoretical Chemistry, University of California, Berkeley, USA.

For assistance, contact: [nikhilagrawal0165@gmail.com](mailto:nikhilagrawal0165@gmail.com).

