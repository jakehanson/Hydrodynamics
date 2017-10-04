![alt text](https://github.com/jakehanson/ISING/blob/dev/Ising_Freeze.gif)

## Overview
This is an Ising Model with a Metpropolis update algorithm. To run the code, you first specify the size of the square lattice (up to 8x8) and the temperature. The model then randomly initializes and evolves the model through a specified number of steps. The lattice size, temeperature, and number of steps can all be modified in **main.cpp**.

There are three primary outputs from the code:
* **params.txt** - The input parameters from main.cpp are output to a file called params.txt.

* **states.txt** - Each spin configuration is assigned a unique integer value specifying its state. The value assigned to a given spin configuration is the integer value of the binary representation of the lattice (e.g. [0 0 1 1] -> 0+0+4+8 = 12). The states.txt file is, therefore, a 1d array specifying the trajectory of the model through the integer valued state space.

* **time_series.txt** - To visualize the evolution of the 2d model requires the entire spin matrix at each point in time. This is what time_series.txt contains. Note, for long simulations of large lattices, the size of this output file can become very large! For example, the maximum lattice size is 8x8, which means a single time step contains 8x8*sizeof(int) = 2^7 Bytes. Therefore, if you run the simulation for a few million time steps, the file size will be roughly 1 GB! To turn off the time_series output, set "spin_output = False" in **main.cpp**.

## Download and Run
To run the code, you first must download the source files **main.cpp**,**functions.cpp*, and **header.h**. Then, you must compile these into an executable:

g++ -std=c++11 main.cpp functions.cpp -o run_sim.exe

Once you have the executable, it can be run from the command line:

./run_sim.exe

**NOTE: EVERY TIME YOU CHANGE A PARAMETER IN MAIN.CPP YOU MUST RECOMPILE THE CODE!!**