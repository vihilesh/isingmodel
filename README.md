Test readme
Check commits
check readme again
To Run IsingModel.Py ( just the ones with PRNG) the command is:
Python IsingModel.py PCG64 --MCSteps 1024 --EQSteps 256 --Lattice 6 --TempPts 16
The command to run Merseene Twister is
Python IsingModel.py MT --MCSteps 1024 --EQSteps 256 --Lattice 6 --TempPts 16

The data should be written to c:\\data\\  - please make a directory before running. 

To Run the quantum model 
python IsingModelQRNG.py Quantum --MCSteps 1024 --EQSteps 256 --Lattice 6 --TempPts 16

The plots and the graph file are stored in c:\\data\\