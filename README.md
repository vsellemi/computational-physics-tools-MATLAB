# Computational Physics in MATLAB
Computational tools for problems in statistical thermodynamics and quantum mechanics, implemented in MATLAB

## Polynomials/Root Finding 

'regulafalsi.m' : implements Regula-Falsi Method to estimate roots of a polynomial

'secant.m' : implements Secant Method to estimate roots of a polynomial

'bisection.m' : implements Bisection Method to estimate roots of a polynomial

'lagrangeinterp.m' : returns the Lagrange interpolating polynomial for function f

## Numerical Solutions of Differential Equations

'numerov.m' implements Numerov's method to solve second order ODEs

'laplace_2d_solutions.m' numerical solution of 2D Laplace equation with various boundary conditions

'diffusion_CN.m' : numerical solution of diffusion equation using Crank-Nicholson algorithm

'numerical_black_scholes.m' : Crank-Nicholson Algorithm to solve the Black-Scholes equation

## Examples

'dielectric_bandstructure.m' : simulated the bandstructure of an infinite stack of dielectric slabs using Metropolis algorithm and finite differences

'hydrogen_wavefunction.m': solves for the ground-state wavefunction and eigenenergy for the hydrogen atom using variational method with Gaussian basis 

'metropolis_hastings_particles.m' : Metropolis-Hastings simulations of diffusion of particles

'qm_scattering.m' : Quantum mechanical scattering simulation using finite differences and transfer matrix method

'quantum_well_solutions.m' : Quantum well numerical eigensolutions

'verlet_anharmonic_osc.m' : simulate anharmonic oscillator using Verlet algorithm

'verlet_damped_pendulum.m' : damped and driven pendulum with Verlet algorithm

'wavefunction_simulations.m' : simulate orbital bandstructure of graphene with finite differences
