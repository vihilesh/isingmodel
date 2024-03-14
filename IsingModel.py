

import argparse
import csv
import io
import os
import pprint
import pandas as pd
import numpy as np
from numpy.random import rand, Generator, PCG64, MT19937, SeedSequence
import matplotlib.pyplot as plt
from scipy.sparse import spdiags,linalg,eye
from enum import Enum



class PRNG(Enum):
    PCG64 = 1
    MT19937 = 2
    QRNG = 3

class RNG :
    def __init__(self, RNG):
        if RNG==PRNG.PCG64 :
            self.prng=Generator(PCG64())
        elif RNG ==  PRNG.MT19937 :
            self.prng = Generator(MT19937())
            

            

def initLattice(RNG, N):

    lattice = RNG.prng.integers(size=(N,N), low=0, high=1, endpoint=True)
    #lattice = 2*np.random.randint(2, size=(N,N))-1
    lattice[lattice==0] =- 1

    return lattice  
    

def mcmove(lattice, RNG, beta, sizeof_lattice):
    '''
    Monte Carlo move using Metropolis algorithm 
    '''
    N = sizeof_lattice
    for i in range(N):
        for j in range(N):
                a =  RNG.prng.integers(low=0, high=N, size=1, dtype=np.int8, endpoint=False)
                b =  RNG.prng.integers(low=0, high=N, size=1, dtype=np.int8, endpoint=False)
                s =  lattice[a, b]
                nb = lattice[(a+1)%N,b] + lattice[a,(b+1)%N] + lattice[(a-1)%N,b] + lattice[a,(b-1)%N]
                cost = 2*s*nb
                
                if cost < 0:
                    s *= -1
                elif  RNG.prng.random() < np.exp(-cost*beta):
                    s *= -1
                lattice[a, b] = s
    return lattice


def calcEnergy(config, sizeof_lattice):
    '''Energy of a given configuration'''
    energy = 0
    N = sizeof_lattice
    for i in range(len(config)):
        for j in range(len(config)):
            S = config[i,j]
            nb = config[(i+1)%N, j] + config[i,(j+1)%N] + config[(i-1)%N, j] + config[i,(j-1)%N]
            energy += -nb*S
    return energy/4

def calcMag(config):
    '''
    Magnetization of a given configuration
    '''
    mag = np.sum(config)
    return mag


def plotfigure(titlestr, E, T, C, X, M):
    f = plt.figure(figsize=(18, 10)); #  
    f.suptitle(titlestr)
    filename = os.path.join( "C:", "\\data", titlestr)

    df = pd.DataFrame({"Temperature": T, "Energy": E, "Magnetisation": M, "Specific Heat":C, "Susceptibility": X})
    #file = str(eqSteps) + "-" + str(mcSteps) + "-" + str(N)
    df.to_csv(filename + ".csv", index=False)

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

    #plt.show()
    #plotstr = rngstr + "-MCS-"+ str(mcSteps) + "-EQS-" + str(eqSteps) + "-Grid-" + str(N) + "-TempPts-" + nt + ".svg"
    plt.savefig(filename+".svg", format="svg")
    plt.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('RNG', choices=["PCG64", "MT"], type=str, help='Random Number Generator')
    parser.add_argument('--MCSteps', type=int, help='Monte Carlo Sweeps')
    parser.add_argument('--EQSteps', type=int, help='Equalibriation steps')
    parser.add_argument('--Lattice', type=int, help='Lattice matrix')
    parser.add_argument('--TempPts', type=int, help='Temperature Points')
    parser.add_argument('--Trial', type=int, help='Trial #')
    args = parser.parse_args()

    print(args.RNG)
    print(args.MCSteps)
    print(args.EQSteps)
    print(args.Lattice)
    print(args.TempPts)
    print(args.Trial)

    nt = args.TempPts
    N = args.Lattice
    mcSteps = args.MCSteps
    eqSteps = args.EQSteps
    rngstr = args.RNG

    if (args.RNG=="PCG64"):
        rng = RNG(PRNG.PCG64)
    elif (args.RNG=="MT") :
        rng = RNG(PRNG.MT19937)
    else :
        raise Exception("Invalid RNG specified")



    trialstr = "-Trial-"+ str(args.Trial) # + "-"
    titlestr = "RNG-" + rngstr + "-MCSteps-" + str(mcSteps) + "-EqSteps-" + str(eqSteps) + "-Ltce-" + str(N) + "-TPts-" + str(args.TempPts)
    titlestr = titlestr+trialstr
    print(titlestr)

    #get a array of floating temp values between the low and high
    T       = np.linspace(1.53, 3.28, nt); 

    #zero out the values for E, M, C and X for nt size. 
    E,M,C,X = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)


    n1, n2  = 1.0/(mcSteps*N*N), 1.0/(mcSteps*mcSteps*N*N) 

    #loop over all temperatues in the nt array
    for tt in range(nt):
        lattice = initLattice(rng, N)        # initialise

        E1 = M1 = E2 = M2 = 0
        iT=1.0/T[tt]; iT2=iT*iT;         # beta in other words
        
        for i in range(eqSteps):         # equilibrate
            mcmove(lattice, rng, iT,N)           # Monte Carlo moves

        for i in range(mcSteps):
            mcmove(lattice, rng, iT,N)           
            Ene = calcEnergy(lattice,N)     # calculate the energy
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
        print (tt)

    plotfigure(titlestr, E, T, C, X, M)

if __name__=='__main__':
    main()

