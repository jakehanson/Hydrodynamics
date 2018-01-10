## SLOPE LIMITERS - STEP FUNCTION
Preserving steep gradients in numerical hydrodynamics is notoriously difficult due to the fact that discrete grids have limited resolution. For example, consider initializing a step function propogating rightward with a fixed velocity. After the first time step, the CFL condition guarantees that the front edge of the step function only PARTIALLY enters the cell on the right. Thus, the new value of this cell will be an average of the step value and 0 (with the exact value depending on how far into the cell the step propogated). This means that we now have a transition cell on the boundary of the step function that is between the step value and 0... or in other words, the clean step has been smeared out. This, in a nutshell, is numerical diffusion.

To combat numerical diffusion, different advection algorithms have been invented. Without going into the details, this directory contains a comparison of some of the most well-known of these algorithms. Note, we evolve a step function through 300 timesteps of dt=0.1 on a grid of 100 cells.

Additional Resources:
http://www2.mpia-hd.mpg.de/~dullemon/lectures/fluiddynamics08/

![alt text](https://github.com/jakehanson/Hydrodynamics/blob/master/SLOPE_LIMITERS/sim.gif)