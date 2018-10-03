export ISS

mutable struct ISS{T<:AbstractFloat} <: Environment
  worldAABBmin::Vector{T}
  worldAABBmax::Vector{T}
  keepin_zones::Vector
  keepout_zones::Vector
  obstacle_set::Vector
end
function ISS{T}() where T
  keepin_zones = Vector{HyperRectangle}(0)
  vars = matread(joinpath(ENV["GUSTO"], "environment","iss.mat"))
  for zone in vars["keepin_zones"]
    push!(keepin_zones,
      HyperRectangle(Vec3f0(zone["corner1"]),Vec3f0(zone["corner2"]-zone["corner1"])))
  end

  worldAABBmin = Inf*ones(T,3)
  worldAABBmax = -Inf*ones(T,3)
  for zone in keepin_zones
    corner1,corner2 = zone.origin, zone.origin+zone.widths
    zone_min = [min(corner1[i], corner2[i]) for i in 1:3]    
    worldAABBmin[worldAABBmin .> zone_min] = zone_min[worldAABBmin .> zone_min]
    zone_max = [max(corner1[i], corner2[i]) for i in 1:3]
    worldAABBmax[worldAABBmax .< zone_max] = zone_max[worldAABBmax .< zone_max]
  end

  keepout_zones = Vector{GeometryTypes.GeometryPrimitive}(0)
  obstacle_set = Vector{GeometryTypes.GeometryPrimitive}(0)

  return ISS{T}(worldAABBmin, worldAABBmax, keepin_zones, keepout_zones, obstacle_set)
end
ISS(::Type{T} = Float64; kwargs...) where {T} = ISS{T}(; kwargs...)

function update_aabb!(env::ISS)
  for zone in env.keepin_zones
    corner1,corner2 = zone.origin, zone.origin+zone.widths
    zone_min = [min(corner1[i], corner2[i]) for i in 1:3]    
    env.worldAABBmin[env.worldAABBmin .> zone_min] = zone_min[env.worldAABBmin .> zone_min]
    zone_max = [max(corner1[i], corner2[i]) for i in 1:3]
    env.worldAABBmax[env.worldAABBmax .< zone_max] = zone_max[env.worldAABBmax .< zone_max]
  end
  warn("Overriding collision world in env()")
end
