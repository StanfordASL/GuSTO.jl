mutable struct Forest{T<:AbstractFloat} <: Environment
  worldAABBmin::Vector{T}
  worldAABBmax::Vector{T}
  keepin_zones::Vector
  keepout_zones::Vector
  obstacle_set::Vector
  obstacles_set_render::Vector
end

function Forest{T}() where T
  vars = matread(joinpath(Pkg.dir("GuSTO"), "src", "environment","forest.mat"))
  
  keepout_zones = Vector{HyperRectangle}(0)
  push!(keepout_zones, HyperRectangle(Vec3f0(-20.,-20.,28.), Vec3f0(160.,160.,20)))
  push!(keepout_zones, HyperRectangle(Vec3f0(-20.,-20.,-20.), Vec3f0(160.,160.,20)))
  push!(keepout_zones, HyperRectangle(Vec3f0(-20.,-20.,-20.), Vec3f0(20.,160.,68)))
  push!(keepout_zones, HyperRectangle(Vec3f0(120.,-20.,-20.), Vec3f0(20.,160.,68)))
  push!(keepout_zones, HyperRectangle(Vec3f0(-20.,-20.,-20.), Vec3f0(160.,20.,68)))
  push!(keepout_zones, HyperRectangle(Vec3f0(-20.,120.,-20.), Vec3f0(160.,20.,68)))

  obstacle_set = Vector{HyperRectangle}(0)
  for idx in 1:2:size(vars["obstacles_infl"],1)
    c1,c2 = vars["obstacles_infl"][idx,:], vars["obstacles_infl"][idx+1,:] 
    push!(obstacle_set,
      HyperRectangle(Vec3f0(c1),Vec3f0(c2-c1)))
  end

  obstacles_set_render = Vector{HyperRectangle}(0)
  keys = ["tree_obstacles", "tower_obstacles"]
  for key in keys
    for idx in 1:2:size(vars[key],1)
      c1,c2 = vars[key][idx,:], vars[key][idx+1,:] 
      push!(obstacles_set_render,
        HyperRectangle(Vec3f0(c1),Vec3f0(c2-c1)))
    end
  end

  worldAABBmin = Inf*ones(T,3)
  worldAABBmax = -Inf*ones(T,3)
  for zone in keepout_zones
    corner1,corner2 = zone.origin, zone.origin+zone.widths
    zone_min = [min(corner1[i], corner2[i]) for i in 1:3]    
    worldAABBmin[worldAABBmin .> zone_min] = zone_min[worldAABBmin .> zone_min]
    zone_max = [max(corner1[i], corner2[i]) for i in 1:3]
    worldAABBmax[worldAABBmax .< zone_max] = zone_max[worldAABBmax .< zone_max]
  end
  
  keepin_zones = Vector{HyperRectangle}(0)
  widths = worldAABBmax-worldAABBmin
  margin = 0.1*widths
  push!(keepin_zones, HyperRectangle(Vec3f0(worldAABBmin-margin),Vec3f0(widths+margin)))

  return Forest{T}(worldAABBmin, worldAABBmax, keepin_zones,
    keepout_zones, obstacle_set, obstacles_set_render)
end
Forest(::Type{T} = Float64; kwargs...) where {T} = Forest{T}(; kwargs...)

# function plot_forest(f::Forest)
#   function plot_shape(rect::HyperRectangle,plot_color::String="blue")
#     lenX,lenY,lenZ = rect.widths 
#   
#     (x1,y1,z1) = rect.origin
#     (x2,y2,z2) = rect.origin+rect.widths
#     
#     xs = [x1; x2; x2; x1; x1]
#     ys = [y1; y1; y2; y2; y1]
#     zs = [z1; z1; z1; z1; z1]
#    
#     zs2 = [z2; z2; z2; z2; z2]
#   
#     xs3 = [x1; x1; x2; x2; x2; x2; x1; x1]
#     ys3 = [y2; y2; y2; y2; y1; y1; y1; y1]
#     zs3 = [z1; z2; z2; z1; z1; z2; z2; z1]
#   
#     PyPlot.mesh(xs, ys, zs, color=plot_color)
#     PyPlot.mesh(xs, ys, zs2, color=plot_color)
#     PyPlot.mesh(xs3, ys3, zs3, color=plot_color)
#   end
# 
#   PyPlot.figure()
#   plot_shape(f.keepin_zones[1],"green")
# 
#   for z in f.keepout_zones
#     plot_shape(z, "red")
#   end
#   for z in f.obstacles_set_render
#     plot_shape(z, "black")
#   end
# end
