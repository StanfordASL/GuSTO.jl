export Table

mutable struct Table{T<:AbstractFloat} <: Environment
  worldAABBmin::Vector{T}
  worldAABBmax::Vector{T}
  keepin_zones::Vector
  keepout_zones::Vector
  obstacle_set::Vector
end

function Table{T}(room::Symbol=:ames) where T
  # 2D granite table

  if room == :ames
    # NASA Ames table (from Andrew Symington)
    # CoM coordinates + robot radius
    radius = 0.15*sqrt(2)
    worldAABBmin, worldAABBmax = [-0.5,-0.75,0.] .- radius, [0.75,0.75,0.001] .+ radius
  else room == :stanford
    # Durand 010 table
    ft2m = 0.3048
    worldAABBmin, worldAABBmax = [0.,0.,0.], [12.,9.,0.001]*ft2m
  end
 
  return Table(worldAABBmin, worldAABBmax)
end
Table(room::Symbol=:ames,::Type{T} = Float64; kwargs...) where {T} = Table{T}(room; kwargs...)

function Table(worldAABBmin::AbstractArray{T}, worldAABBmax::AbstractArray{T}) where T
  keepin_zones = Vector{HyperRectangle}(undef, 0)
  (length(worldAABBmin) == 2) ? push!(worldAABBmin, T(-1)) : nothing
  (length(worldAABBmax) == 2) ? push!(worldAABBmax, T(1)) : nothing
  push!(keepin_zones, HyperRectangle(worldAABBmin..., (worldAABBmax-worldAABBmin)...))

  # Add obstacles surrounding table
  keepout_zones = Vector{GeometryTypes.GeometryPrimitive}(undef, 0)
  koz_min = []
  koz_max = []
  a = 10.  # keep-out zone size (should be larger than 1.05)
  push!(koz_min, [worldAABBmax[1],   -a, -a])  # x
  push!(koz_max, [worldAABBmax[1]+a,  a,  a])
  push!(koz_min, [worldAABBmin[1]-a, -a, -a])  # -x
  push!(koz_max, [worldAABBmin[1],    a,  a])
  push!(koz_min, [-a, worldAABBmax[2],   -a])  # y
  push!(koz_max, [ a, worldAABBmax[2]+a,  a])
  push!(koz_min, [-a, worldAABBmin[2]-a, -a])  # -y
  push!(koz_max, [ a, worldAABBmin[2],    a])

  for i in 1:length(koz_min)
    push!(keepout_zones, HyperRectangle(koz_min[i]..., (koz_max[i]-koz_min[i])...))
  end
  obstacle_set = Vector{GeometryTypes.GeometryPrimitive}(undef, 0)

  return Table{T}(worldAABBmin, worldAABBmax, keepin_zones, keepout_zones, obstacle_set)
end
