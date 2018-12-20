export Astrobee2D

mutable struct Astrobee2D{T<:AbstractFloat} <: Robot
  mass::T
  J::T
  Jinv::T
  n_thrusters::Int
  s::T
  r::T
  hard_limit_vel::T
  hard_limit_accel::T 
  hard_limit_ω::T
  hard_limit_α::T
  btCollisionObject

  xb::Vector{T}
  Jcollision
end
function Astrobee2D{T}() where T
  n_thrusters = 12
  s = 0.5*0.305   # each side of cube is 30.5cm
  r = sqrt(2)*s   # inflate to sphere

  # ground robot param: freeflyer/astrobee/config/worlds/granite.config
  mass = 14.4
  hard_limit_vel = 0.20 
  hard_limit_accel = 0.02
  hard_limit_ω = 10*π/180 
  hard_limit_α = 10*π/180 
 
  J = 0.1083
  Jinv = inv(J)

  xb = [0.; 0.15; 0.]
  Jcollision = []

  # btCollisionObject = BulletCollision.convex_hull_cylinder(SVector{3}(zeros(3)), SVector{3}(0.,0.,1.), r)
  btCollisionObjects = BulletCollision.BulletCollisionObjectPtr[]
  push!(btCollisionObjects, BulletCollision.convex_hull_cylinder(SVector{3}(zeros(3)), SVector{3}(0.,0.,1.), r))
  push!(btCollisionObjects, BulletCollision.convex_hull_cylinder(SVector{3}(xb...), SVector{3}(0.,0.,1.), r))

  btCollisionObject = BulletCollision.compound_collision_object(btCollisionObjects)

  # new astrobee instance
  return Astrobee2D{T}(mass, J, Jinv, n_thrusters, s, r, hard_limit_vel, hard_limit_accel, hard_limit_ω, hard_limit_α, btCollisionObject, xb, Jcollision)
end
Astrobee2D(::Type{T} = Float64; kwargs...) where {T} = Astrobee2D{T}(; kwargs...)
