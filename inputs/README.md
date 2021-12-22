# Inputs

Modify the file `inputs.txt` to change the parameters used in the calculation. The parameters names are written in 
the left hand column, and the parameter values are in the right hand column, separated by spaces. At the moment,
only J (the coupling constant) can take multiple values. All other parameters must have single numerical values
or the booleans "True" or "False". Lines beginning with a `#` character are ignored.

The module `util.py` contains utility functions that are used to interpret the contents of the `inputs.txt` file.

## Parameter descriptions

Parameter | Description
--------- | ----------
N | Grid size (square 2D N x N lattice)
T | Temperatures in Kelvin, written (lowest temperature):(highest temperature):(temperature interval)
e | Numerical value of electronic charge, in SI units
J | Coupling constant for interaction term of Ising Hamiltonian, in Joules
Nit | The number of iterations of the Metropolis algorithm to be used in the calculation at each temperature
error | True or False, depending on whether or not to calculate error bars. The error calculation is done by computing the total magnetisation of the system (expensive) every "oft" calculations (see next table entry), storing the last "n_last" results (see metropolis function in `computations.py`) and calculating their standard deviation.
oft | The number of iterations in between each of the calculations of the total energy and magnetisation of the system.