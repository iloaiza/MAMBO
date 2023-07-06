# THIS REPOSITORY HAS BEEN DEPRECATED. GO TO [`QuantumMAMBO.jl`](https://github.com/iloaiza/QuantumMAMBO.jl) FOR THE UP-TO-DATE VERSION.

# MAMBO: Efficient many-body routines in Julia

## Using MAMBO
MAMBO includes efficient implementations for fermionic and qubit operators. To obtain results shown in Ref.(1), run on a terminal (e.g. for LiH):

julia L1.jl lih

All options and tolerances can be seen in config.jl. Make sure that the python environment which was used to build PyCall on Julia is active when running this command (see install.sh script for more info).

## Installation
All required packages and installation info can be seen/installed in install.sh.

Fast installation: execute install.sh in a terminal. Set PY_INSTALL=true in file to install local python environment with necessary packages. Can use custom directory for julia packages by uncommenting JL_DIR lines (useful for computing environments where writing to folder with packages is not allowed while running calculations, e.g. Niagara on Compute Canada)
Requires installing julia packages (can be done by accessing the julia package manager in a julia session with ']', then adding packages listed on install.sh. Can also be installed by running the "install.sh" script on a terminal, set PY_INSTALL to true(false) to do(not) install python environment along with julia packages. Make sure to build PyCall with correct python environment (check install.sh script, or PyCall github page for more info).

Requires/creates a python executable with installed packages:
'pip install pyscf openfermion openfermionpyscf h5py'

The python virtual environment should then be activated before running MAMBO routines with
'source VIRTUAL_ENVIRONMENT_DIRECTORY/bin/activate'
and making sure that julia is running with the correct JULIA_DEPOT_PATH bash variable if not using default package installation directory.

## Module overview

### UTILS folder
	- bliss.jl: functions for BLISS routine (see Ref. 2)
	- cost.jl: functions for calculating different norms of operators. Mainly 1- and 2-norms of fermionic operators.
	- decompose.jl: CSA, DF, and related decompositions of fermionic operators
	- ferm_utils.py: utilities for fermionic operators in python, interfaces with openfermion
	- fermionic.jl: utilities for MAMBO fermionic operators class (F_OP, defined in structures.jl)
	- guesses.jl: initial guesses for decomposition routines
	- ham_utils.py: python utilities for the electronic structure Hamiltonian, interfaces with openfermion
	- lcu.jl: calculation of lcu 1-norms for different decompositions
	- linprog.jl: linear programming routines for symmetry-shift optimization (see Refs. 1 and 2, corresponds to "partial" routine in Ref. 2)
	- majorana.jl: utilities for MAMBO Majorana operators class (M_OP, defined in structures.jl)
	- orbitals.jl: orbital optimization routine for 1-norm minimization (see Koridon et al., Phys. Rev. Res. 3 (3), 2021. Material is also covered in Refs. 1 and 2)
	- parallel.jl: code with parallel capabilities, mainly for trotter bounds (full implementation is under progress)
	- projectors.jl: builds projectors of Fock space into constant number of electrons subspace, useful for Trotter bounds (under progress)
	- py_qubits.jl: python utilities for qubit operators, interfaces with openfermion
	- py_utils.jl: julia interface to all python modules and openfermion
	- qubit.jl: utilities for MAMBO qubit operators class (Q_OP, defined in structures.jl)
	- saving.jl: save-load utilities for decompositions and optimization results, uses HDF5
	- structures.jl: definition of classes for many-body operators
	- symmetries.jl: building of symmetry operators, e.g. Sz, Ne
	- symplectic.jl: utilities for representing and manipulating qubit space as symplectic vectors
	- trotter.jl: Trotterization implementation, errors and bounds (under progress)
	- unitaries.jl: unitary transformations related to fermionic MAMBO operators (F_OP)
	- wrappers.jl: runner functions which run workflows for obtaining all necessary quantities for e.g. tables in Refs. 1 and 2

### Main folder
	- L1.jl: full workflow for obtaining LCU 1-norms for all decompositions/methods
	- build.jl: single module for loading all utilities
	- config.jl: general parameters and settings for all functions
	- install.sh: installer file, see Installation section for more information


## References
This code was developped and used for all results in the publications:

[1] - I. Loaiza, A. Marefat Khah, N. Wiebe, and A. F. Izmaylov, Reducing molecular electronic Hamiltonian simulation cost for Linear Combination of Unitaries approaches. arXiv:2208.08272 (2022).

[2] - I. Loaiza, A. F. Izmaylov, Reducing the molecular electronic Hamiltonian encoding costs on quantum computers by symmetry shifts. arXiv:2304.13772 (2023).
