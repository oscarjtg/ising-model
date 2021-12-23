# Module for Ising Model computations

import numpy as np

debug = True
debug2 = False
debug3 = False


def grid_generator(N):
    """
    Function to generate an NxN grid of 1s and (-1)s
    representing spin up and spin down respectively.
    """
    grid = np.zeros((N, N))
    #print(grid)

    # Generate random numbers, uniform distribution in [0,1].
    random_numbers = np.random.uniform(size=N**2) 
    if debug2: print(random_numbers)
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

    if debug2: print(grid)
    if debug2: print('sum of elements in grid is ', count)
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


def metropolis_iteration(grid, i, j, N, J, beta, random_number):
    """ Carries out steps 2-5 of Metropolis algorithm described in README. """
    Ei = nearest_neighbour_interaction(grid, i, j, N, J)
    if debug2: print('before', grid[i,j]); print('Ei = ', Ei)
    
    grid[i, j] = -1.0 * grid[i, j]
    Ef = nearest_neighbour_interaction(grid, i, j, N, J)
    dE = Ef - Ei
    if debug2: 
        print('after', grid[i, j])
        print('Ef = ', Ef)
        print('dE = ', dE)
    
    if dE < 0: # Flipping spin leads to reduction in energy, i.e. is favourable
        if debug3: print('dE<0, remain flipped')
        gridf = grid
        if debug3: print('initial', gridf[i, j]); print('final', gridf[i, j])
    else:
        B_factor = np.exp(-beta * dE)
        if debug3:
            print('beta =', beta)
            print('B-factor=', B_factor)
            print('random number is', random_number)
        if B_factor > random_number:
            if debug3: print('dE>0 but remains flipped')
            gridf = grid
            if debug3: print('initial', gridf[i, j]); print('final', gridf[i, j])
        else:
            if debug3: print('dE > 0, flip back')
            grid[i, j] = -1.0*grid[i, j]
            gridf = grid
            if debug3: print('initial', gridf[i, j]); print('final', gridf[i, j])

    if debug2 or debug3: print('------')
    return gridf


def metropolis(sim, grid1, beta, progplot):
    """ 
    Generates random numbers to select lattice site locations and calls
    metropolis_iteration function.
    """
    N = sim.N
    Nit = sim.Nit
    J = sim.J
    HN = sim.HN
    error = sim.error
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
        grid1 = metropolis_iteration(grid1, u_rand1[k], u_rand2[k], N, J, 
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


def run_calculations(sim):
    """
    Reads Simulation class object "sim".
    For each temperature in sim.temperatures:
    (1) Calculates beta factor (1/kB T) at the temperature
    (2) Calls grid_generator(N) to generate NxN 2D lattice of randomly-distributed up and down spins.
    (3) Calls metropolis() to carry out the Metropolis algorithm, returning magnetisation and uncertainty
    after sim.Nit iteration steps.
    
    Returns arrays of magnetisation and magnetisation uncertainty (arbitrary units) 
    containing the simulation results at each temperature.
    """
    number_temperatures = len(sim.temperatures)
    magnetisation = np.linspace(0, 0, number_temperatures)
    magnetisation_error = np.linspace(0, 0, number_temperatures)
    for temperature in sim.temperatures:
        beta = 1/(sim.kB * temperature)
        grid = grid_generator(sim.N)
        n = np.where(sim.temperatures == temperature)
        magnetisation[n], magnetisation_error[n] = metropolis(sim, grid, beta, False)

    return magnetisation, magnetisation_error
