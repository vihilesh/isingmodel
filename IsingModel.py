

import numpy as np
from numpy.random import rand, Generator, PCG64, MT19937, SeedSequence
import matplotlib.pyplot as plt
from scipy.sparse import spdiags,linalg,eye
from enum import Enum


nt      = 16          #  number of temperature points
N       = 6           #  size of the lattice, N x N
eqSteps = 2**8        #  number of MC sweeps for equilibration
mcSteps = 2**10       #  number of MC sweeps for calculation



class PRNG(Enum):
    PCG64 = 1
    MT19937 = 2

class RNG :
    def __init__(self, RNG):
        if RNG==PRNG.PCG64 :
            self.prng=Generator(PCG64())
        elif RNG ==  PRNG.MT19937 :
            self.prng = Generator(MT19937())

def initLattice(RNG, N):

    lattice = RNG.prng.integers(size=(N,N), low=0, high=1, endpoint=True)
    lattice[lattice==0] =- 1

    return lattice  
    


def mcmove(lattice, RNG, beta):
    '''
    Monte Carlo move using Metropolis algorithm 
    '''
    for i in range(N):
        for j in range(N):
                a =  RNG.prng.integers(low=0, high=N, size=1, dtype=np.int8, endpoint=False)
                b =  RNG.prng.integers(low=0, high=N, size=1, dtype=np.int8, endpoint=False)
                a = a[0]
                b = b[0]
                s =  lattice[a, b]
                nb = lattice[(a+1)%N,b] + lattice[a,(b+1)%N] + lattice[(a-1)%N,b] + lattice[a,(b-1)%N]
                cost = 2*s*nb
                
                if cost < 0:
                    s *= -1
                elif  RNG.prng.random() < np.exp(-cost*beta):
                    s *= -1
                lattice[a, b] = s
    return lattice

def calcEnergy(config):
    '''
    Energy of a given configuration
    '''
    energy = 0 
    
    for i in range(len(config)):
        for j in range(len(config)):
            S = config[i,j]
            nb = config[(i+1)%N, j] + config[i,(j+1)%N] + config[(i-1)%N, j] + config[i,(j-1)%N]
            energy += -nb*S
    return energy/2.  # to compensate for over-counting

def calcMag(config):
    '''
    Magnetization of a given configuration
    '''
    mag = np.sum(config)
    return mag


#rng = RNG(PRNG.PCG64)
rng = RNG(PRNG.MT19937)



T       = np.linspace(1.53, 3.28, nt); 
E,M,C,X = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)
n1, n2  = 1.0/(mcSteps*N*N), 1.0/(mcSteps*mcSteps*N*N) 

for tt in range(nt):
    lattice = initLattice(rng, N)        # initialise

    E1 = M1 = E2 = M2 = 0
    iT=1.0/T[tt]; iT2=iT*iT;
    
    for i in range(eqSteps):         # equilibrate
        mcmove(lattice, rng, iT)           # Monte Carlo moves

    for i in range(mcSteps):
        mcmove(lattice, rng, iT)           
        Ene = calcEnergy(lattice)     # calculate the energy
        Mag = calcMag(lattice)        # calculate the magnetisation

        E1 = E1 + Ene
        M1 = M1 + Mag
        M2 = M2 + Mag*Mag 
        E2 = E2 + Ene*Ene


    # divide by number of sites and iteractions to obtain intensive values    
    E[tt] = n1*E1
    M[tt] = n1*M1
    C[tt] = (n1*E2 - n2*E1*E1)*iT2
    X[tt] = (n1*M2 - n2*M1*M1)*iT

    f = plt.figure(figsize=(18, 10)); #  



sp =  f.add_subplot(2, 2, 1 );
plt.scatter(T, E, s=50, marker='o', color='IndianRed')
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Energy ", fontsize=20);         plt.axis('tight');


sp =  f.add_subplot(2, 2, 2 );
plt.scatter(T, abs(M), s=50, marker='o', color='RoyalBlue')
plt.xlabel("Temperature (T)", fontsize=20); 
plt.ylabel("Magnetization ", fontsize=20);   plt.axis('tight');


sp =  f.add_subplot(2, 2, 3 );
plt.scatter(T, C, s=50, marker='o', color='IndianRed')
plt.xlabel("Temperature (T)", fontsize=20);  
plt.ylabel("Specific Heat ", fontsize=20);   plt.axis('tight');   


sp =  f.add_subplot(2, 2, 4 );
plt.scatter(T, X, s=50, marker='o', color='RoyalBlue')
plt.xlabel("Temperature (T)", fontsize=20); 
plt.ylabel("Susceptibility", fontsize=20);   plt.axis('tight');

plt.show()
