# GuSTO.jl: Guaranteed Sequential Trajectory Optimization

This is a Julia suite for trajectory optimization using the Guaranteed Sequential Trajectory Optimization (GuSTO) framework. Details can be found in [this paper](https://arxiv.org/abs/1903.00155).

GuSTO.jl runs on julia v1.0+, though an older version running on [julia v0.6.4](https://julialang.org/downloads/oldreleases.html) can be found in the julia-v0.6 branch.

Also required are the [BulletCollision.jl](https://github.com/StanfordASL/BulletCollision.jl) and [AstrobeeRobot.jl](https://github.com/StanfordASL/AstrobeeRobot.jl) packages. GuSTO.jl performs optimization through the [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) interface, and Gurobi and Ipopt are currently used in examples.

# Quickstart
An example notebook can be run through:
```
jupyter notebook examples/freeflyerSE2.ipynb 
```

Click to watch demo video:

<a href="https://www.youtube.com/watch?v=GHehE-If5nY" target="_blank"><img src="https://img.youtube.com/vi/GHehE-If5nY/0.jpg" 
alt="cla" width="240" height="180" border="10" /></a>
