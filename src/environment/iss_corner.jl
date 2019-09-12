export ISSCorner

mutable struct ISSCorner{T<:AbstractFloat} <: Environment
  worldAABBmin::Vector{T}
  worldAABBmax::Vector{T}
  keepin_zones::Vector
  keepout_zones::Vector
  obstacle_set::Vector
end

function ISSCorner{T}() where T
  keepin_zones = Vector{HyperRectangle}(0)
  vars = matread(joinpath(Pkg.dir("GuSTO"), "src", "environment","iss_corner.mat"))
  for zone in vars["keepin_zones"]
    push!(keepin_zones,
      HyperRectangle(Vec3f0(zone["corner1"][:]),Vec3f0(zone["corner2"][:]-zone["corner1"][:])))
  end

  keepout_zones = Vector{GeometryTypes.GeometryPrimitive}(0)
  for zone in vars["keepout_zones"]
    push!(keepout_zones,
      HyperRectangle(Vec3f0(zone["corner1"][:]),Vec3f0(zone["corner2"][:]-zone["corner1"][:])))
  end

  obstacle_set = Vector{GeometryTypes.GeometryPrimitive}(0)

  worldAABBmin = Inf*ones(T,3)
  worldAABBmax = -Inf*ones(T,3)
  for zone in (keepin_zones..., keepout_zones...)
    corner1, corner2 = zone.origin, zone.origin + zone.widths
    zone_min = [min(corner1[i], corner2[i]) for i in 1:3]    
    worldAABBmin[worldAABBmin .> zone_min] = zone_min[worldAABBmin .> zone_min]
    zone_max = [max(corner1[i], corner2[i]) for i in 1:3]
    worldAABBmax[worldAABBmax .< zone_max] = zone_max[worldAABBmax .< zone_max]
  end

  return ISSCorner{T}(worldAABBmin, worldAABBmax, keepin_zones, keepout_zones, obstacle_set)
end
ISSCorner(::Type{T} = Float64; kwargs...) where {T} = ISSCorner{T}(; kwargs...)

function update_aabb!(env::ISSCorner)
  for zone in env.keepin_zones
    corner1, corner2 = zone.origin, zone.origin + zone.widths
    zone_min = [min(corner1[i], corner2[i]) for i in 1:3]    
    env.worldAABBmin[env.worldAABBmin .> zone_min] = zone_min[env.worldAABBmin .> zone_min]
    zone_max = [max(corner1[i], corner2[i]) for i in 1:3]
    env.worldAABBmax[env.worldAABBmax .< zone_max] = zone_max[env.worldAABBmax .< zone_max]
  end
  warn("Overriding collision world in env()")
end

function add_obstacles1!(env::ISSCorner)
  vars = matread(joinpath(Pkg.dir("GuSTO"), "src", "environment","iss_corner.mat"))
  for zone in vars["rectangles"]
    push!(env.obstacle_set,
      HyperRectangle(Vec3f0(zone["corner1"][:]),Vec3f0(zone["corner2"][:]-zone["corner1"][:])))
  end

  for zone in vars["spheres"]
    push!(env.obstacle_set,
      HyperSphere(Point3f0(zone["center"][:]), Float32(zone["radius"])))
  end
end



function add_obstacles2!(env::ISSCorner)
  # ----------
  # RECTANGLES

  bottom_left = Vec3f0(+8.50,-0.50,+4.5)
  top_right   = Vec3f0(+9.25,+0.25,+5.25)
  push!(env.obstacle_set, HyperRectangle(bottom_left, top_right-bottom_left))
  bottom_left = Vec3f0(+10.50,-0.5,+4.75)
  top_right   = Vec3f0(+11.00,-0.0,+5.00)
  push!(env.obstacle_set, HyperRectangle(bottom_left, top_right-bottom_left))
  bottom_left = Vec3f0(+10.75,+3.5,+4.50)
  top_right   = Vec3f0(+11.50,+4.2,+5.25)
  push!(env.obstacle_set, HyperRectangle(bottom_left, top_right-bottom_left))
  bottom_left = Vec3f0(+10.00,+5.5,+4.75)
  top_right   = Vec3f0(+10.75,+7.0,+5.75)
  push!(env.obstacle_set, HyperRectangle(bottom_left, top_right-bottom_left))
  # ----------


  # vars = matread(joinpath(Pkg.dir("GuSTO"), "src", "environment","iss_corner.mat"))
  # for zone in vars["rectangles"]
  #   push!(env.obstacle_set,
  #     HyperRectangle(Vec3f0(zone["corner1"][:]),Vec3f0(zone["corner2"][:]-zone["corner1"][:])))
  # end

  # for zone in vars["spheres"]
  #   push!(env.obstacle_set,
  #     HyperSphere(Point3f0(zone["center"][:]), Float32(zone["radius"])))
  # end
end

function add_obstacles3!(env::ISSCorner)
  # ----------
  # RECTANGLES

  bottom_left = Vec3f0(+8.50,-0.50,+4.5)
  top_right   = Vec3f0(+9.25,+0.25,+5.25)
  push!(env.obstacle_set, HyperRectangle(bottom_left, top_right-bottom_left))
  bottom_left = Vec3f0(+10.50,-0.5,+4.75)
  top_right   = Vec3f0(+11.00,-0.0,+5.00)
  push!(env.obstacle_set, HyperRectangle(bottom_left, top_right-bottom_left))
  bottom_left = [+10.75 ; +3.5 ; +4.50 ]
  top_right   = [+11.50 ; +4.2 ; +5.25 ]
  center = ( bottom_left + top_right ) ./ 2.
  radius = 0.4
  push!(env.obstacle_set, HyperSphere(Point3f0(center), Float32(radius)))
  bottom_left = Vec3f0(+10.00,+5.5,+4.75)
  top_right   = Vec3f0(+10.75,+7.0,+5.75)
  push!(env.obstacle_set, HyperRectangle(bottom_left, top_right-bottom_left))
  # ----------


  # vars = matread(joinpath(Pkg.dir("GuSTO"), "src", "environment","iss_corner.mat"))
  # for zone in vars["rectangles"]
  #   push!(env.obstacle_set,
  #     HyperRectangle(Vec3f0(zone["corner1"][:]),Vec3f0(zone["corner2"][:]-zone["corner1"][:])))
  # end

  # for zone in vars["spheres"]
  #   push!(env.obstacle_set,
  #     HyperSphere(Point3f0(zone["center"][:]), Float32(zone["radius"])))
  # end
end
