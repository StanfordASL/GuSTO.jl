# GuSTO.jl: Guaranteed Sequential Trajectory Optimization

This is a Julia suite for trajectory optimization using the Guaranteed Sequential Trajectory Optimization (GuSTO) framework. Details can be found in [this paper](https://arxiv.org/abs/1903.00155).

GuSTO.jl runs on julia v1.X, although an older version running on [julia v0.6.4](https://julialang.org/downloads/oldreleases.html) can be found in the julia-v0.6 branch.

Also required are the [BulletCollision.jl](https://github.com/StanfordASL/BulletCollision.jl) and [AstrobeeRobot.jl](https://github.com/StanfordASL/AstrobeeRobot.jl) packages. GuSTO.jl performs optimization through the [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) interface, and Gurobi and Ipopt are currently used in examples.

## Quickstart
An example notebook can be run through:
```
jupyter notebook examples/freeflyerSE2.ipynb 
```

Click to watch demo video:

<a href="https://www.youtube.com/watch?v=GHehE-If5nY" target="_blank"><img src="https://img.youtube.com/vi/GHehE-If5nY/0.jpg" 
alt="cla" width="240" height="180" border="10" /></a>

## References
* R. Bonalli, A. Cauligi, A. Bylard, and M. Pavone, ["GuSTO: Guaranteed Sequential Trajectory Optimization via Sequential Convex Programming,"](https://arxiv.org/abs/1903.00155) in *Proc. IEEE Conf. on Robotics and Automation*, 2019.
* R. Bonalli, A. Bylard, A. Cauligi, T. Lew, and M. Pavone, ["Trajectory Optimization on Manifolds: A Theoretically-Guaranteed Embedded Sequential Convex Programming Approach,"](https://arxiv.org/abs/1903.00155) in *Robotics: Science and Systems*, 2019.
* R. Bonalli, T. Lew, and M. Pavone, ["Analysis of Theoretical and Numerical Properties of Sequential Convex Programming for Continuous-Time Optimal Control,"](https://arxiv.org/abs/2009.05038) in *IEEE Transactions on Automatic Control*, 2022.

## Update
We recommend using the Julia implementation of GuSTO available at [https://github.com/UW-ACL/SCPToolbox.jl](https://github.com/UW-ACL/SCPToolbox.jl).
