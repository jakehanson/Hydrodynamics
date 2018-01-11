## ALGORITHM COMPARISON - STEP FUNCTIONS AND NUMERICAL DIFFUSION
Preserving steep gradients in numerical hydrodynamics is notoriously difficult due to the fact that discrete grids have limited resolution. For example, consider initializing a step function propogating rightward with a fixed velocity. After the first time step, the CFL condition guarantees that the front edge of the step function only PARTIALLY enters the cell on the right. Thus, the new value of this cell will be an average of the step value and 0 (with the exact value depending on how far into the cell the step propogated). This means that we now have a transition cell on the boundary of the step function that is between the somewhere between the step value and 0... or in other words, the clean step has been smeared out. This, in a nutshell, is numerical diffusion.

To combat numerical diffusion, different advection algorithms have been invented. Without going into the details, this directory contains a comparison of some of the most well-known of these algorithms. We look at slope limiters, which are higher order algorithms designed to prevent overshoot of slope estimation, and fllux limiters, which have result in smooth second order flux-conserved advection on smooth parts of the solution and flux-conserving first order donor-cell advection near sharp or steep gradients. For shock-capturing hydrodynamics, flux limiters are essential, although they have the potential to steepen existing gradients in a non-physical way.

Additional Resources:
http://www2.mpia-hd.mpg.de/~dullemon/lectures/fluiddynamics08/

In the image, we evolve a step function through 300 timesteps of dt=0.1 on a grid of 100 cells using various slope limiters.

![alt text](https://github.com/jakehanson/Hydrodynamics/blob/master/ALGORITHM_COMPARISON/SLOPE_LIMITERS/sim.gif)