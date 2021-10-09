# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 11:01:28 2021

@author: Oscar Tovey Garcia

Ising model on a 10x10 grid using the Metropolis algorithm.
"""
import sys

import numpy as np
import matplotlib.pyplot as plt
import time

debug = True


def read_inputs(filename):
    """Reads the input file and creates a dict with appropriate entries."""
    input_parameters = {}
    with open(filename, 'r') as f:
        file_text = f.read()
    file_lines = file_text.split('\n')
    for line in file_lines:
        if debug: print(line)
        if line.startswith('#') or line.startswith('\n'):
            if debug: print(f"skipped line {line}")
            continue

        pair = line.split()
        if debug: print(pair)
        if len(pair) < 2:
            if pair == []:
                if debug: print("expect empty list: {pair}")
                continue
            raise ValueError(f"read_input: too few values in line - {pair}")
        elif len(pair) > 2:
            raise ValueError(f"read_inputs: too many values in line - {pair}")
        
        key = pair[0]
        value = pair[1]

        try:
            input_parameters[key] = float(value)
        except ValueError:
            if value == "True":
                input_parameters[key] = True
            elif value == "False":
                input_parameters[key] = False
            else:
                input_parameters[key] = value

    for key in input_parameters.keys():
        value = input_parameters[key]
        if isinstance(value, str):
            if value.endswith('{e}'):
                if debug: print(value)
                try:
                    electronic_charge = input_parameters['e']
                except KeyError:
                    raise KeyError(f"Have not specified electronic charge in {filename}")
                factor = float(value[:-3])
                if debug:
                    print(f'electronic charge = {electronic_charge}')
                    print(f'factor = {factor}')

                input_parameters[key] = factor * electronic_charge

    if debug:
        print(input_parameters)
    return input_parameters



def grid_generator(N):
    """
    Function to generate an NxN grid of 1s and (-1)s
    representing spin up and spin down respectively.
    """
    grid = np.zeros((N, N))
    #print(grid)

    # Generate random numbers, uniform distribution in [0,1].
    random_numbers = np.random.uniform(size=N**2) 
    #print(random_numbers) # Uncomment to print the random numbers.
    count = 0

    # This for loop randomly allocates +1 or -1 to a grid entry
    for i in range(N):
        for j in range(N):
            if random_numbers[i + j] > 0.5:
                grid[i, j] = 1
                count += 1
            else:
                grid[i, j] = -1
                count -= 1

    # print(grid)
    # print('sum of elements in grid is ', c)
    return grid


def nearest_neighbour_interaction(grid, i, j, N, J=1):
    """
    Function to calculate interaction of nearest neighbours to spin (i, j) in 
    grid with periodic boundary conditions (so e.g. if grid is 10x10, 
    element (0, 0) couples to (1, 0), (9, 0), (0, 1) and (0, 9)).
    """
    grid_point = grid[i, j]
    #print(int((i-1)%N))
    #print(i)
    # Nearest neighbours:
    nn1 = grid[int((i - 1) % N), j]
    nn2 = grid[int((i + 1) % N), j]
    nn3 = grid[i, int((j - 1) % N)]
    nn4 = grid[i, int((j + 1) % N)]
    nearest_neighbour_sum = nn1 + nn2 + nn3 + nn4

    return -J * grid_point * nearest_neighbour_sum


def full_interaction_energy(grid, N=10, J=1):
    """
    Calculates the Hamiltonian with coupling constant J.
    """
    H_sum = 0 # Initialise Hamiltonian sum.
    
    # Sum all the nearest neighbour interaction energies.
    for i in range(N):
        for j in range(N):
            H_sum += nearest_neighbour_interaction(grid, i, j, N, J)
            
    H = 0.5 * H_sum # Factor of 0.5 to remove double counting.
    return H


def Metropolis_iteration(grid, i, j, J, B, ur):
    #print('before', grid[i,j])
    Ei = nearest_neighbour_interaction(grid, i, j, N, J)
    #print('Ei = ', Ei)
    grid[i, j] = -1.0 * grid[i, j]
    #print('after', grid[i, j])
    Ef = nearest_neighbour_interaction(grid, i, j, N, J)
    #print('Ef = ', Ef)
    dE = Ef - Ei
    #print('dE = ', dE)
    if dE < 0: # Flipping spin leads to reduction in energy, i.e. is favourable
        #print('dE<0, remain flipped')
        gridf = grid
        #print(gridf[i, j])
        #print('final', gridf[i, j])
    else:
        #ur = np.random.uniform(1)
        B_factor = np.exp(-B * dE)
        #print('beta =',B)
        #print('B-factor=',B_factor)
        #print('random number is',ur)
        if B_factor > ur:
            #print('dE>0 but remains flipped')
            gridf = grid
            #print('final',gridf[i,j])
        else:
            #print('dE > 0, flip back')
            grid[i, j] = -1.0*grid[i, j]
            gridf = grid
            #print('final', gridf[i, j])

    #print('------')
    return gridf


def Metropolis(grid1, N, Nit, J, beta, progplot, HN, error):
    u_rand1 = np.random.randint(0, N, size=Nit)
    u_rand2 = np.random.randint(0, N, size=Nit)
    u_rand3 = np.random.uniform(size=Nit)
    n_last = 1000;

    if error is True:
        # Initialise array storing the last magnetisation values
        # for subsequent calculation of standard deviation.
        if Nit > n_last:
            #H_last = np.zeros(n_last)
            M_last = np.zeros(n_last)
        else:
            #H_last = np.zeros(Nit)
            M_last = np.zeros(Nit)

    #print(u_rand1)
    #print(u_rand2)
    #print(u_rand3)

    k = 0
    l = 0
    m = 0
    
    if progplot is True:
        H = np.linspace(HN + 1)
    
    for k in range(Nit):
        grid1 = Metropolis_iteration(grid1, u_rand1[k], u_rand2[k], J, 
                                     beta, u_rand3[k])
        if k % (Nit / HN) == 0 and progplot is True:
            H[l + 1] = full_interaction_energy(grid1, N, J)
            l += 1
        # store the last 200 values to calculate mean and standard dev (for errorbar plot)
        if error is True and (Nit-k) <= n_last:
            #H_last[m] = full_interaction_energy(grid1,N,J)
            M_last[m] = np.sum(grid1)
            m += 1
        k += 1

    if progplot is True:
        print('final grid', grid1)

        plt.figure(20)
        t = np.linspace(0,20,HN+1)
        plt.plot(t,H)

    if error is True:
        #Mag = np.mean(H_last)
        #Mag_dev = np.std(H_last)
        Mag = np.mean(M_last)
        Mag_dev = np.std(M_last)
    else:
        #Mag = full_interaction_energy(grid1,N,J)
        Mag = np.sum(grid1)
        Mag_dev = 0
    #print(grid1)
    return Mag, Mag_dev


def calculate_magnetisation(temperatures, J):
    number_temperatures = len(temperatures)
    magnetisation = np.linspace(0, 0, number_temperatures)
    magnetisation_error = np.linspace(0, 0, number_temperatures)
    for temperature in temperatures:
        grid = grid_generator(N)
        beta = 1/(kB * temperature)
        n = np.where(temperatures == temperature)
        magnetisation[n], magnetisation_error[n] = Metropolis(grid, N, Nit, J, beta, False, HN, error)

    return magnetisation, magnetisation_error



if __name__ == "__main__":
    if len(sys.argv) > 2:
        filename = sys.argv[2]
    else:
        filename = 'inputs.txt'

    params = read_inputs(filename)

    N = int(params['N'])
    #T = params['T']
    kB = params['kB']
    #beta = 1 / (kB * T) # Thermodynamic beta.
    e = params['e']
    J = params['J']
    Nit = int(params['Nit'])
    oft = int(params['oft'])
    HN = int(Nit / oft)
    error = params['error']
    
    temperatures = np.arange(100, 900, 50)
    #temperatures = [0.001, 0.01, 0.1, 1, 10, 30, 60, 80, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 5000, 8000, 10000, 20000, 50000, 1.0e5]
    number_temperatures = len(temperatures)

    if debug:
        print(f"temperatures = {temperatures}")
        print(f"number of temperatures = {number_temperatures}")

    J1 = J
    tic1 = time.process_time()
    magnetisation_1, magnetisation_1_error = calculate_magnetisation(temperatures, J1)
    toc1 = time.process_time()
    if debug: print(f'Calculation 1 time is {toc1 - tic1} s')

    J2 = 2 * J
    tic2 = time.process_time()
    magnetisation_2, magnetisation_2_error = calculate_magnetisation(temperatures, J2)
    toc2 = time.process_time()
    if debug: print(f'Calculation 2 time is {toc2 - tic2} s')

    fig, ax = plt.subplots()
    ax.errorbar(temperatures, np.abs(magnetisation_1), magnetisation_1_error, 0, 
        'bs', ms=4, label=f'J = {J1} meV')
    ax.errorbar(temperatures, np.abs(magnetisation_2), magnetisation_2_error, 0, 
        'mo', ms=4, label=f'J = {J2} meV')

    plt.xlabel('Temperature / K')
    plt.ylabel('|Magnetisation| / au')
    plt.title('Ising Model')
    ax.legend(loc='upper right')

    plt.savefig('Ising_model5.png')
    plt.savefig('Ising_model5.pdf')

    toc5 = time.process_time()
    if debug: print(f'Total runtime is {toc5 - tic1} s.')

    plt.show()