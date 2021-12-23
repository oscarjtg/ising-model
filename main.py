# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 11:01:28 2021

@author: Oscar Tovey Garcia

Ising model on a 10x10 grid using the metropolis algorithm.
"""
import sys

import numpy as np
import time

import computations as comps
import inputs.util
import plots

debug = True
debug2 = False
debug3 = False

class Simulation:
    """ 
    Class that runs simulations of Ising Model using Metropolis Algorithm,
    storing both inputs and results. 
    """
    def __init__(self, inputfile):
        params = inputs.util.read_inputs(filename)

        self.N = int(params['N'])
        self.temperatures = params['T']
        self.kB = params['kB']
        self.e = params['e']
        self.Nit = int(params['Nit'])
        self.oft = int(params['oft'])
        self.HN = int(self.Nit / self.oft)
        self.error = params['error']
        self.number_temperatures = len(self.temperatures)

    def set_J(self, J):
        self.J = J

    def magnetisation_calculation(self):
        mag, err = comps.run_calculations(self)
        self.magnetisation = mag
        self.magnetisation_error = err


if __name__ == "__main__":
    
    # Process any command-line arguments.
    if debug2: print(sys.argv)
    if len(sys.argv) > 1:
        filename = sys.argv[1]
        if debug: print('Used command line input')
    else:
        filename = 'inputs/inputs.txt'
    if debug: print(filename)

    # Read the input file and extract parameters into a dict called "params".
    params = inputs.util.read_inputs(filename, debug=True)
    
    # Process the coupling constant input(s) into a Python list.
    J = params['J']
    if isinstance(J, str):
        J_list = list(float(J))
    elif isinstance(J, list):
        J_list = J
    else:
        raise ValueError(f'J from params is neither a str nor list: it is {J}.')
    if debug: print(J_list)
    
    # Carry out simulation for each J in the input file, timing each calculation.
    results = {}
    for J in J_list:
        num = J_list.index(J) + 1
        tic = time.process_time()
        results[num] = Simulation(filename)
        results[num].set_J(J)
        results[num].magnetisation_calculation()
        toc = time.process_time()
        if debug: print(f'Calculation {num} time is {toc - tic} s')
    
    # Generate one plot of magnetisation vs temperature, with different points of each J.
    plots.magnetisation_temperature(results, debug)

    