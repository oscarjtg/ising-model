# Ising Model

Uses the Metropolis algorithm to numerically solve the Ising model for the net magnetisation of a ferromagnet at different temperatures.

A plot of |magnetisation| against temperature is produced, comparing systems with different coupling constants J. 
The ferromagnet-to-paramagnet phase transition is observed at a higher temperature when the coupling constant is larger.

## Instructions

To run the code, first check `inputs.txt` in the inputs directory and modify the parameter values as desired. 
A description of the parameters is found in the README.md file in the inputs directory.

Then navigate back to the ising-model directory and input ```python ising_model.py``` at the command line.
The computation will proceed, and a plot will pop out as a separate file. This is automatically saved into 
the output directory as a PDF and as a PNG (`ising_model.pdf` and `ising_model.png`).

## Further description

The module called `computations.py` contains the source code (in Python) that produces a 2D spin lattice, iterates
using the Metropolis algorithm and calculates the magnetisation in the final state. The algorithm works as follows:
1. Choose a lattice site at random.
2. Calculate the nearest-neighbour interaction energy, Ei, from the interaction term of the Ising Hamiltonian using nearest-neighbour spins of the chosen lattice site only.
3. Flip the spin of the chosen lattice site and repeat the nearest neighbour interaction energy calculation, to get Ef.
4. If Ef < Ei, keep the spin flipped and proceed to step 6.
5. If Ef > Ei, calculate the associated Boltzmann factor exp(-(Ef-Ei)/Kb T) and compare it to a random number in between 0 and 1. If the random number is larger than the Boltzmann factor, flip the spin back to how it was previously and proceed to step 6. Otherwise, keep the spin flipped and proceed to step 6.
6. Repeat the iteration from step 1.

## Potential improvements

Using the current parameters, the program tests at three different coupling constants, each requiring roughly
5 seconds to compute. The total run time of the script is less than 20 seconds. However, it does not always 
finish converging in this time, so longer iteration time would be required for more accurate results.

The main improvements are:

1. Modify the code to stop computation for each temperature only once the system's magnetisation has equilibrated, rather than simply after a fixed number of iterations.
2. Produce an animation showing the evolution of the system as the algorithm reaches equilibrium at a given temperature.