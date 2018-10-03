#__precompile__()
# __precompile__(false)
# module GuSTO

using BulletCollision

using PyPlot
using RigidBodySim, RigidBodyDynamics
using MeshCat, MeshCatMechanisms
using CoordinateTransformations
using Interact, Reactive
using MAT
using FileIO
using MeshIO
using MechanismGeometries
using ForwardDiff
import GeometryTypes: HyperRectangle, HyperSphere, HomogenousMesh, Vec, Vec3f0, Point, Point3f0
import ColorTypes: RGB, RGBA

using StaticArrays
using StaticArrays.FixedSizeArrays
using JuMP, Ipopt, MathProgBase
using Convex, Gurobi, Mosek, SCS

using NLsolve, DifferentialEquations

using MAT
using GeometryTypes

include("types.jl")

include("traj_opt.jl")
include("scp.jl")
include("shooting.jl")

include("utils/quat_functions.jl")

include("robot.jl")
include("dynamics.jl")
include("environment.jl")

# end
