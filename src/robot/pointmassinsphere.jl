export PointMassInSphere

# Simple double integrator
mutable struct PointMassInSphere{T<:AbstractFloat} <: Robot
  mass::T
  r::T
  hard_limit_vel::T
  hard_limit_accel::T 
  btCollisionObject
end
function PointMassInSphere{T}() where T
  # Default values
  mass = 1.0
  r = 0.1
  hard_limit_vel = 1.0
  hard_limit_accel = 1.0

  btCollisionObject = BulletCollision.convex_hull_cylinder(SVector{3}(zeros(3)), SVector{3}(0.,0.,1.), r)

  return PointMassInSphere{T}(mass, r, hard_limit_vel, hard_limit_accel, btCollisionObject)
end
PointMassInSphere(::Type{T} = Float64; kwargs...) where {T} = PointMassInSphere{T}(; kwargs...)
