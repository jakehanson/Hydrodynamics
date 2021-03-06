## ALGORITHM COMPARISON - STEP FUNCTIONS AND NUMERICAL DIFFUSION
Preserving steep gradients in numerical hydrodynamics is notoriously difficult due to the fact that discrete grids have limited resolution. For example, consider initializing a step function propogating rightward with a fixed velocity. After the first time step, the CFL condition guarantees that the front edge of the step function only PARTIALLY enters the cell on the right. Thus, the new value of this cell will be an average of the step value and 0 (with the exact value depending on how far into the cell the step propagated). This means that we now have a transition cell on the boundary of the step function that is somewhere between the step value and 0... or in other words, the clean step has been smeared out. This, in a nutshell, is numerical diffusion.

To combat numerical diffusion, different advection algorithms have been invented. Without going into the details, this directory contains a comparison of some of the more well-known of these algorithms. We look at slope limiters, which are higher order algorithms designed to prevent overshoot of slope estimation, and flux limiters, which result in smooth second order flux-conserved advection on smooth parts of the solution and flux-conserving first order donor-cell advection near sharp or steep gradients. For shock-capturing hydrodynamics, flux limiters are essential, although they have the potential to steepen existing gradients in a non-physical way.

Additional Resources:
http://www2.mpia-hd.mpg.de/~dullemon/lectures/fluiddynamics08/

In the first image, we evolve a step function through 300 timesteps of dt=0.1 on a grid of 100 cells using various slope limiters. We see that without slope estimation (Donor-Cell) the step function suffers strong numerical diffusion, while downwind slope estimation (Lax-Wendroff) overshoots and becomes numerically unstable near steep gradients. Thus, a slope limiter (Superbee) limits numerical diffusion while ensuring that we don't overshoot.

The second image implements various flux limiters. We see that they steepen gradients and further limit numerical diffusion.

![alt text](https://github.com/jakehanson/Hydrodynamics/blob/master/ALGORITHM_COMPARISON/SLOPE_LIMITERS/sim.gif)
![alt text](https://github.com/jakehanson/Hydrodynamics/blob/master/ALGORITHM_COMPARISON/FLUX_LIMITERS/sim.gif)
