# __precompile__(false)

using BulletCollision

using Plots
using RigidBodySim, RigidBodyDynamics
using MeshCat, MeshCatMechanisms
using CoordinateTransformations
using Interact, Reactive
using MAT
using FileIO
using MeshIO
using MechanismGeometries
using ForwardDiff
import GeometryTypes: HyperRectangle, HyperSphere, HomogenousMesh
import ColorTypes: RGB, RGBA

using StaticArrays
using StaticArrays.FixedSizeArrays
using JuMP, Convex
using Ipopt, Mosek, SCS
using MathOptInterfaceMosek
using Gurobi
using MathProgBase, MathOptInterface
using NLsolve, DifferentialEquations

using MAT
using GeometryTypes
using FillArrays
using LinearAlgebra

const MOI = MathOptInterface

include("types.jl")

include("traj_opt.jl")
include("scp.jl")
include("shooting.jl")

include("utils/quat_functions.jl")

include("robot.jl")
include("dynamics.jl")
include("environment.jl")

