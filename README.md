# IntroHeatEquation.jl
Workshop on methods to simulate and control simple thermal dynamics.


## Content

This repository contains [Pluto.jl](https://github.com/fonsp/Pluto.jl) scripts to simulate the free and controlled one-dimensional heat equation with various boundary conditions.

### Modeling
The free (uncontrolled) heat equation models with various boundary condidtions (Dirichlet, Neumann, Robin) are presented here.

- Dirichlet: constant temperature on both boundaries
  - [analytical](https://github.com/stephans3/IntroHeatEquation.jl/blob/main/src/modeling/dirichlet_analytical.jl) and 
  - [numerical](https://github.com/stephans3/IntroHeatEquation.jl/blob/main/src/modeling/dirichlet_numerical.jl) solution
- zero Neumann: insulated sides / no flux towards environment
  - [analytical](https://github.com/stephans3/IntroHeatEquation.jl/blob/main/src/modeling/neumann_analytical.jl) and 
  - [numerical](https://github.com/stephans3/IntroHeatEquation.jl/blob/main/src/modeling/neumann_numerical.jl) solution
- Nonlinear: heat transfer and heat radiation towards environment 
  - [numerical](https://github.com/stephans3/IntroHeatEquation.jl/blob/main/src/modeling/robin_numerical.jl) solution

### Control
Simple control approaches for the discretized heat equation are presented here. 

- [Proportional control](https://github.com/stephans3/IntroHeatEquation.jl/blob/main/src/control/prop_control.jl)
- [PI control](https://github.com/stephans3/IntroHeatEquation.jl/blob/main/src/control/pi_control.jl) (advanced)


## How to use it
1. Install [Julia](https://julialang.org/)
2. Add Pluto.jl to your Julia environment
3. Start a Pluto session
4. Insert the link of simulation scripts from this repository in the Pluto session, see the textbox "Open from file" 


Please note: This is work in progress.
