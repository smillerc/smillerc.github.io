# About

<!-- ![](/assets/images/profile_pic.png) -->

<!-- ![](__site/assets/images/profile_pic.png) -->


I'm a computational physicist working at the [Laboratory for Laser Energetics](https://www.lle.rochester.edu/) at the University of Rochester, NY. I did my Ph.D. in Mechanical Engineering specializing in hydrodynamic instabilities that occur in inertial confinement fusion. More specifically, my research included high-fidelity simulations of micron-scale internal defects in cryogenic targets and described how these defects evolve and create seeds for hydrodynamic instability growth. In my day-to-day job I maintain and develop high performance scientific code that simulate laser direct-drive inertial confinement fusion implosions in 1D, 2D, and 3D. 

## Interests

- Modern code development practices
- Julia, modern Fortran, C/C++
- High performance CFD methods

## Projects

 - [MPIHaloArrays.jl](https://github.com/smillerc/MPIHaloArrays.jl) -- A new `Array` type in Julia that manages the halo-exchange process (via MPI) that commonly occurs in large-scale domain-decomposed simulation codes
 - [Cato](https://github.com/smillerc/cato) -- A flexible modern Fortran (2018+) code that solves the Euler fluid equations within the finite volume method. Parallelization is accomplished using Coarrays and OpenMP.
 - [coarray_field](https://github.com/smillerc/coarray_field) -- A simple high-level field object in Fortran 2018 that manages domain decomposition via coarrays