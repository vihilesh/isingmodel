

import argparse
import csv
import io
import linecache
import pprint
import pandas as pd
import numpy as np
from numpy.random import rand, Generator, PCG64, MT19937, SeedSequence
import matplotlib.pyplot as plt
from scipy.sparse import spdiags,linalg,eye
from enum import Enum


nt      = 16          #  number of temperature points
N       = 6           #  size of the lattice, N x N
eqSteps = 2**8        #  number of MC sweeps for equilibration
mcSteps = 2**10     #  number of MC sweeps for calculation



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
    lattice[lattice==0] =- 1

    return lattice  

def get_Integer_QRNG_from_file(QRNG_int_file, QRNG_int_file_num_lines):
    line1 = int(linecache.getline("QRNGIntegers.txt", np.random.randint(1, QRNG_int_file_num_lines)).strip())
    line2 = int(linecache.getline("QRNGIntegers.txt", np.random.randint(1,QRNG_int_file_num_lines)).strip())
    return line1, line2

def get_Float_QRNG_from_file(QRNG_float_file, QRNG_float_file_num_lines):
    try:
        line_in_file = np.random.randint(1,QRNG_float_file_num_lines)
        line_text = linecache.getline("QRNGFloat.txt", line_in_file).strip()
        line = float(line_text)
        return line
    except:
        print(line_in_file)
    


def mcmove(lattice, RNG, beta, QRNG_int_file, QRNG_int_file_num_lines, QRNG_float_file, QRNG_num_lines):
    '''
    Monte Carlo move using Metropolis algorithm 
    '''
    for i in range(N):
        for j in range(N):
                a, b =  get_Integer_QRNG_from_file(QRNG_int_file, QRNG_int_file_num_lines)
                s =  lattice[a, b]
                nb = lattice[(a+1)%N,b] + lattice[a,(b+1)%N] + lattice[(a-1)%N,b] + lattice[a,(b-1)%N]
                cost = 2*s*nb
                
                if cost < 0:
                    s *= -1
#                elif  get_Float_QRNG_from_file(QRNG_float_file, QRNG_num_lines) < np.exp(-cost*beta):
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
rngstr = ""
# if (PRNG.MT19937):
#     rng = RNG(PRNG.MT19937)
#     rngstr = "Mersenne Twister"

if (PRNG.PCG64):
    rng = RNG(PRNG.PCG64)
    rngstr = "PCG64"

titlestr = "RNG-" + "QRNG" + "-MCSteps-" + str(mcSteps) + "-EqSteps-" + str(eqSteps) + "-Lattice:" + str(N)


parser = argparse.ArgumentParser(prog="IsingModel")
subparser = parser.add_subparsers(title="actions")
parser_list = subparser.add_parser('RNG choices', help="Specify a RNG")
parser_list.add_argument("-RNG", required=True, choices=['PCG64','MT19937'])
args = parser.parse_args()
print(vars(args))


#get a array of floating temp values between the low and high
T       = np.linspace(1.53, 3.28, nt); 

#zero out the values for E, M, C and X for nt size. 
E,M,C,X = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)


n1, n2  = 1.0/(mcSteps*N*N), 1.0/(mcSteps*mcSteps*N*N) 

QRNG_int_file_linecount = 0
QRNG_Int_File =  open("QRNGIntegers.txt", 'rb')
QRNG_int_file_linecount = sum(1 for _ in QRNG_Int_File)

QRNG_float_file_linecount = 0
QRNG_float_file =  open("QRNGFloat.txt", 'rb')
QRNG_float_file_linecount = sum(1 for _ in QRNG_float_file)

QRNG_Int_File.close
QRNG_float_file.close

#loop over all temperatues in the nt array
for tt in range(nt):
    lattice = initLattice(rng, N)        # initialise

    E1 = M1 = E2 = M2 = 0
    iT=1.0/T[tt]; iT2=iT*iT;         # beta in other words
    
    for i in range(eqSteps):         # equilibrate
        mcmove(lattice, rng, iT, QRNG_Int_File, QRNG_int_file_linecount, QRNG_float_file,QRNG_float_file_linecount)           # Monte Carlo moves

    for i in range(mcSteps):
        mcmove(lattice, rng, iT, QRNG_Int_File, QRNG_int_file_linecount, QRNG_float_file, QRNG_float_file_linecount)           
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
f.suptitle(titlestr)

df = pd.DataFrame({"Temperature": T, "Energy": E, "Magnetisation": M, "Specific Heat":C, "Susceptibility": X})
file = str(eqSteps) + "-" + str(mcSteps) + "-" + str(N)
df.to_csv(file+".csv", index=False)

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
plotstr = rngstr + "-MCS-"+ str(mcSteps) + "-EQS-" + str(eqSteps) + "-Grid-" + str(N) + ".svg"
plt.savefig(plotstr, format="svg")
plt.close()

QRNG_Int_File.close
QRNG_float_file.close

