# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 11:01:28 2021

@author: OscarTG
"""

# Ising model on a 10x10 grid using the Metropolis algorithm yay!

import numpy as np
import matplotlib.pyplot as plt
import time

# declaring parameters
N = 10 # grid size
T = 100 # temperature in Kelvin
kB = 1.381e-23 # Boltzmann constant in m^2 kg s^-2 K^-1
#kB = 1
beta = 1/(kB * T) # thermodynamic beta
e = 1.602e-19 # electronic charge
#e = 1
J = 0.01*e # coupling constant in J aka kg m^2 s^-2
Nit = 100000 # number of iterations
oft = 500 # how often to calculate H
HN = int(Nit/oft)
error = 'y' # plot error bars y/n

# function to generate an NxN grid of 1s and (-1)s
# representing spin up and spin down respecitvely
def grid_generator(N):
    grid = np.zeros((N,N))
    #print(grid)

    a = np.random.uniform(size=N**2) #generate random numbers, uniform distribution in [0,1]
    #print(a) #uncomment to print the random numbers
    c = 0

    # this for loop randomly allocates +1 or -1 to a grid entry
    for i in range(N):
        for j in range(N):
            if a[i+j] > 0.5: #compare random number to 0.5 to decide if array element should be +1 or -1
                grid[i,j] = 1
                c += 1
            else:
                grid[i,j] = -1
                c -= 1

    # print(grid)
    # print('sum of elements in grid is ', c) 
    # c tells us the sum of all the elements in the grid, so we can quickly see if there are more 1s or -1s

    return grid


# function to calculate interaction of nearest neighbours to spin (i,j) in grid
# with periodic boundary conditions (so e.g. element (0,0) couple to (1,0), (9,0), (0,1) and (0,9)
def nearest_neighbour_interaction(grid,i,j,J=1):
    Hij = grid[(i-1)%10,j%10] + grid[(i-1)%10,j] + grid[(i+1)%10,j] + grid[i,(j-1)%10] + grid[i,(j+1)%10]
    return -J*grid[i,j]*Hij


# function to calculate the Hamiltonian with coupling constant J 
def full_interaction_energy(grid, N=10, J=1):
    H_sum = 0 #initialise our Hamiltonian sum
    
    # sum all the nearest neighbour interaction energies
    for i in range(N):
        for j in range(N):
            H_sum += nearest_neighbour_interaction(grid,i,j,J)
            
    H = 0.5*H_sum # factor of 0.5 to remove double counting
    return H

def Metropolis_iteration(grid,i,j,J,B,ur):
    #print('before', grid[i,j])
    Ei = nearest_neighbour_interaction(grid,i,j,J)
    #print('Ei = ', Ei)
    grid[i,j] = -1.0*grid[i,j]
    #print('after', grid[i,j])
    Ef = nearest_neighbour_interaction(grid,i,j,J)
    #print('Ef = ', Ef)
    dE = Ef - Ei
    #print('dE = ', dE)
    if dE < 0: # flipping spin leads to reduction in energy, i.e. is favourable
        #print('dE<0, remain flipped')
        gridf = grid
        #print(gridf[i,j])
        #print('final',gridf[i,j])
    else:
        #ur = np.random.uniform(1)
        B_factor = np.exp(-B*dE)
        #print('beta =',B)
        #print('B-factor=',B_factor)
        #print('random number is',ur)
        if B_factor > ur:
            #print('dE>0 but remains flipped')
            gridf = grid
            #print('final',gridf[i,j])
        else:
            #print('dE>0, flip back')
            grid[i,j] = -1.0*grid[i,j]
            gridf = grid
            #print('final',gridf[i,j])

    #print('------')
    return gridf


tic = time.process_time()
#grid1 = grid_generator(N)
#toc = time.process_time()
#print(grid1)
#toc2 = time.process_time()
#print('grid_generator time is ', (toc-tic)*1000, 'ms')
#print('grid_generator time + printing is ', (toc2-tic)*1000, 'ms')


#H = np.linspace(0,0,HN+1)
#H[0] = full_interaction_energy(grid1,N,J)
#print(H)
#progplot = 'n'

def Metropolis(grid1,N,Nit,J,beta,progplot,HN,error):
    u_rand1 = np.random.randint(0,N,size=Nit)
    u_rand2 = np.random.randint(0,N,size=Nit)
    u_rand3 = np.random.uniform(size=Nit)
    n_last = 1000;

    if error == 'y':
        if Nit > n_last:
            #H_last = np.zeros(n_last)
            M_last = np.zeros(n_last)
        else:
            #H_last = np.zeros(Nit)
            M_last = np.zeros(Nit)

    #print(u_rand1)
    #print(u_rand2)
    #print(u_rand3)

    k=0
    l=0
    m=0
    
    if progplot == 'y':
        H = np.linspace(HN+1)
    
    for k in range(Nit):
        grid1 = Metropolis_iteration(grid1,u_rand1[k],u_rand2[k],J,beta,u_rand3[k])
        if k % (Nit/HN) == 0 and progplot == 'y':
            H[l+1] = full_interaction_energy(grid1,N,J)
            l += 1
        # store the last 200 values to calculate mean and standard dev (for errorbar plot)
        if error == 'y' and (Nit-k) <= n_last:
            #H_last[m] = full_interaction_energy(grid1,N,J)
            M_last[m] = np.sum(grid1)
            m += 1
        k += 1

    if progplot == 'y':
        print('final grid', grid1)

        plt.figure(20)
        t = np.linspace(0,20,HN+1)
        plt.plot(t,H)

    if error == 'y':
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


#Metropolis(grid1,N,Nit,J,beta,'y',HN)   

"""
a = np.array([10,30])
b = np.arange(60,200,20)
c = np.arange(220,800,20)
d = np.arange(650,1000,50)
e = np.arange(1100,1500,100)

#1D array imbedding

T = np.zeros(a.size + b.size + c.size + d.size + e.size)
T[:a.size] = a
T[a.size:(a.size+b.size)] = b
T[a.size+b.size:(a.size+b.size+c.size)] = c
T[a.size+b.size+c.size:(a.size+b.size+c.size+d.size)] = d
T[a.size+b.size+c.size+d.size:] = e
"""

T = np.arange(100,1000,100)
#T = [100,200,300,400,500]


#T = [0.001, 0.01, 0.1, 1, 10, 30, 60, 80, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 5000, 8000, 10000, 20000, 50000, 1.0e5]
#T = np.linspace(100,1000,100)
print(T)
M1 = np.linspace(0,0,len(T))
M_dev1 = np.linspace(0,0,len(T))
p=0
for p in range(len(T)):
    gridp = grid_generator(N)
    beta = 1/(kB * T[p]) 
    M1[p], M_dev1[p] = Metropolis(gridp,N,Nit,J,beta,'n',HN,error)
    p += 1


fig, ax = plt.subplots()
ax.errorbar(T,M1,M_dev1,0,'bs',ms=4,label='J = 10 meV')

toc3 = time.process_time()
print('Calculation 1 time is ', (toc3-tic), 's')
tic4 = time.process_time()

# 2nd run
M2 = np.linspace(0,0,len(T))
M_dev2 = np.linspace(0,0,len(T))
p=0
for p in range(len(T)):
    gridp = grid_generator(N)
    beta = 1/(kB * T[p]) 
    M2[p], M_dev2[p] = Metropolis(gridp,N,Nit,J*2,beta,'n',HN,error)
    p += 1

toc4 = time.process_time()
print('Calculation 2 time is ', (toc3-tic), 's')

ax.errorbar(T,M2,M_dev2,0,'mo',ms=4,label='J = 20 meV')

plt.xlabel('Temperature / K')
plt.ylabel('Magnetisation / au')
plt.title('Ising Model')
ax.legend(loc='upper right')

plt.savefig('Ising_model5.png')
plt.savefig('Ising_model5.pdf')

toc5 = time.process_time()
print('Total runtime is', (toc5-tic), 's')

plt.show()








