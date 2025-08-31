

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

## numerical_param.py

This is an input file to specify numerical parameters like the number of grid points, tolerance criteria, mixing ratios for non-linear solvers, ncc cutoffs for _Dedalus_, etc. Note that the equations being solved here are highly non-linear and hence some amount of tuning of these numerical parameters is required for efficient calculations.

## physical_param.py 

This is an input file to specify physical environment variables like salt concentrations, ion valencies, born radii of ions, excluded volumes of ions and solvent, domain size, dielectric permittivity, temperature, etc. In the second part of this file are derived non-dimensional variables from these input parameters, all the calculations are done in these non-dimensional variables. 

## dh_2plate.py

solves the linearized mean-field Poisson-Boltzmann or Debye-Hueckel theory for electrical double layers. The solution of this is usually used as an initial guess to solve for full mean-field PB in pb_2plate.py.

## pb_2plate.py

Solves the full mean-field Poisson-Boltzmann theory for electrical double layers between two charged plates. This is a non-linear boundary value problem whose initial guess can come from dh_2plate.py or another solution for mean-field PB. The solution of this can be used as an initial guess to solve for the modified Gaussian renormalized fluctuation theory in mgrf_2plate.py. 

## mgrf_2plate.py

calculates the solution to modified Gaussian renormalized fluctuation theory for a 2-plate system. This is also a non-linear boundary value problem whose initial guess can come from pb_2plate.py or a solution for mgrf with another set of parameters. This function requires various properties like screening lengths, concentration profiles, self-energies, etc. Functions for these properties are described below.

## num_concn.py
Has three functions. nconc_mgrf calculates the coefficient in front of the exp(-z\psi) in the mgrf_2plate.py. It also outputs the coefficients that are needed to calculate the Jacobian for the Newton-raphson iterative scheme. nconc_complete is the function to calculate the concentration profile for a given psi profile with n_initial as the initial guess. nconc_pb calculates number density profiles for mean-field PB. 

## selfe_2plate.py

This file includes functions to calculate the self-energy profiles for a double layer for a 2-plate system based on the equations given in supplemental material of Agrawal and Wang, Phys. Rev. Lett. 2022, 129, 228001. The functions in this file use another file called greens_function_2plate.py which evaluates the fourier transform of the green's functions in the interface.

## selfe_bulk.py

This file includes functions to calculate the self-energy for the bulk solution based on the equations given in supplemental material of Agrawal and Wang, Phys. Rev. Lett. 2022, 129, 228001. Note that although there is an analytical solution for self-energy in the bulk we calculate it numerically to cancel out the numerical errors between self-energy in interface and the bulk. The functions in this file use another file called greens_function_bulk.py which evaluates the fourier transform of the green's functions in the bulk.

## greens_function_2plate.py and greens_function_bulk.py

File to calculate Fourier transforms of G and Go in the interface and bulk respectively.

## calculate.py

This file contains functions to evaluate properties like screening length, ionic strength, incompressibility fields, and charge density profiles. There is also a function called interpolator to interpolate electrostatic potential and ion density profiles to increase or decrease grid points. A function to calculate the residual of Gauss law is also given.

## energy_2plate.py

Functions to calculate the grand free energy of the interface and bulk for both mean-field PB as well as modified Gaussian renormalized fluctuation theory.

## simulator_pb.py

This code saves the solution of the modified Gaussian renormalized fluctuation theory in a .h5 file for the input parameters given in the two *_param.py files. This file uses the solution of pb_2plate.py as the initial guess to solve for mgrf_2plate.py. The input variables for pb_2plate.py and mgrf_2plate can be set separately using the file physical_param.py.

## simulator.py

This code saves the solution of the modified Gaussian renormalized fluctuation theory in a .h5 file for the input parameters given in the two *_param.py files. This file uses a saved solution of mgrf_2plate.py as the initial guess. The physical variables for this saved solution and the final parameters for which we want the double layer structure can be set separately using the file physical_param.py. The variables deciding which saved solution to choose as initial guess end with "_in_d", for ex: sigma_in_d.

## packages.py

This Python file contains the import statements for all the Python libraries that are needed for this package. We suggest that you create a separate conda environment where all these libraries are installed.

## Running the code

The code can be run using any of the following commands based on your needs: 

python simulator.py physical_param.py

python simulator_pb.py physical_param.py

Note that numerical_param.py has been directly imported into the relavant .py files.

## Contact:
This code was developed by Nikhil Agrawal in the lab of Prof. Rui Wang, Pitzer Center for Theoretical Chemistry, University of California, Berkeley, USA. If you need any help feel free to write to Nikhil at nikhilagrawal0165@gmail.com.  

