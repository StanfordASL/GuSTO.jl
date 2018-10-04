GuSTO.jl
========

This is a suite for trajectory optimization using the Guaranteed Sequential Trajectory Optimization (GuSTO) framework. Details can be found in [this paper](http://asl.stanford.edu/wp-content/papercite-data/pdf/Bonalli.Cauligi.Bylard.Pavone.ICRA19.pdf).

Installation
------------
The code has been tested for Ubuntu 16.04 and requires [julia v0.6.4](https://julialang.org/downloads/oldreleases.html). To install GuSTO, run the following command from the Julia REPL:

```
Pkg.clone("git://github.com/StanfordASL/GuSTO.jl.git")
```

The Julia interfaces to the required optimization solvers will be installed during the build step, but the user must check the specific package instructions for downloading the binary or license files as necessary:
 * The convex optimization problems are designed to be solved with either [Gurobi](https://github.com/JuliaOpt/Gurobi.jl) or [Mosek](https://github.com/JuliaOpt/Mosek.jl), but the user can specify the free [SCS](https://github.com/JuliaOpt/SCS.jl) solver as well. 
 *  The shooting method scripts require [Ipopt](https://github.com/JuliaOpt/Ipopt.jl).

Quickstart
------------
An example notebook can be run through:
```
jupyter notebook examples/astrobeeSE3.ipynb 
```

Click to watch demo video:

<a href="https://www.youtube.com/watch?v=GHehE-If5nY" target="_blank"><img src="https://img.youtube.com/vi/GHehE-If5nY/0.jpg" 
alt="cla" width="240" height="180" border="10" /></a>
