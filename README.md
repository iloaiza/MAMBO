# MAMBO: Efficient many-body routines in Julia

## Using MAMBO
MAMBO includes efficient implementations for fermionic and qubit operators. To obtain results shown in Ref.(1), run on a terminal (e.g. for LiH using 10 parallel processes):

julia -p 10 L1.jl lih

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

## References
This code was developped and used for all results in the publication:

[1] - I. Loaiza, A. Marefat Khah, N. Wiebe, and A. F. Izmaylov, Reducing molecular electronic Hamiltonian simulation cost for Linear Combination of Unitaries approaches. arXiv.2208.08272 (2022).

