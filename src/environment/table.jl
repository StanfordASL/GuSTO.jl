export Table

mutable struct Table{T<:AbstractFloat} <: Environment
  worldAABBmin::Vector{T}
  worldAABBmax::Vector{T}
  keepin_zones::Vector
  keepout_zones::Vector
  obstacle_set::Vector

  # -----------------------------------
  # CoRL- Additional start / end  zones
  safe_region_type::String # "ball" or "rectangle" (2d)
  # for each dim "i", provides a constraint as
  # [coeffs[i][1]*x + coeffs[i][2]*y + coeffs[i][3]*z + coeffs[i][4] < 0
  safe_region_coeffs::Array{Vector{T}} 
  safe_region_smallest_width::T
  safe_region_height::T
  safe_region_width::T

  # goal_region_type::String # "ball" or "rectangle" (2d)
  # goal_region_position::Vector{T}
  # goal_region_radius::T
  # -----------------------------------
end

function Table{T}(room::Symbol=:ames) where T
  # 2D granite table

  if room == :ames
    # NASA Ames table (from Andrew Symington)
    # CoM coordinates + robot radius
    radius = 0.15*sqrt(2)
    worldAABBmin, worldAABBmax = [-0.5,-0.75,0.]-radius, [0.75,0.75,0.001]+radius
  else room == :stanford
    # Durand 010 table
    ft2m = 0.3048
    worldAABBmin, worldAABBmax = [0.,0.,0.], [12.,9.,0.001]*ft2m
  end
 
  return Table(worldAABBmin, worldAABBmax)
end
Table(room::Symbol=:ames,::Type{T} = Float64; kwargs...) where {T} = Table{T}(room; kwargs...)

function Table{T}(worldAABBmin::AbstractArray{T}, worldAABBmax::AbstractArray{T})
  keepin_zones = Vector{HyperRectangle}(0)
  (length(worldAABBmin) == 2) ? push!(worldAABBmin, T(-1)) : nothing
  (length(worldAABBmax) == 2) ? push!(worldAABBmax, T(1)) : nothing
  push!(keepin_zones, HyperRectangle(worldAABBmin..., (worldAABBmax-worldAABBmin)...))

  # Add obstacles surrounding table
  keepout_zones = Vector{GeometryTypes.GeometryPrimitive}(0)
  koz_min = []
  koz_max = []

  B_original = false
  if B_original
    # ------------------
    # Keepin - stay in workspace -> constraints to not come out of it
    a = 10.  # keep-out zone size (should be larger than 1.05)
    push!(koz_min, [worldAABBmax[1],   -a, -a])  # x
    push!(koz_max, [worldAABBmax[1]+a,  a,  a])
    push!(koz_min, [worldAABBmin[1]-a, -a, -a])  # -x
    push!(koz_max, [worldAABBmin[1],    a,  a])
    push!(koz_min, [-a, worldAABBmax[2],   -a])  # y
    push!(koz_max, [ a, worldAABBmax[2]+a,  a])
    push!(koz_min, [-a, worldAABBmin[2]-a, -a])  # -y
    push!(koz_max, [ a, worldAABBmin[2],    a])
    # ------------------

  else
    # ------------------
    # Keepin - stay in workspace -> constraints to not come out of it
    a = 10.  # keep-out zone size (should be larger than 1.05)
    push!(koz_min, [worldAABBmax[1],   -a, -a])  # x
    push!(koz_max, [worldAABBmax[1]+a,  a,  a])
    push!(koz_min, [worldAABBmin[1]-a, -a, -a])  # -x
    push!(koz_max, [worldAABBmin[1],    a,  a])
    push!(koz_min, [-a, worldAABBmax[2],   -a])  # y
    push!(koz_max, [ a, worldAABBmax[2]+a,  a])
    push!(koz_min, [-a, worldAABBmin[2]-a, -a])  # -y
    push!(koz_max, [ a, worldAABBmin[2],    a])


    # --------------------------
    # Canyon (passage in middle)
    middle = worldAABBmin + 0.5*(worldAABBmax-worldAABBmin)
    a,b,c = 0.25, 0.31, 0.01
    a,b,c = 0.25, 0.35, 0.01
    # a,b,c = 0.25, 0.41, 0.01
    # a,b,c = 0.25, 0.51, 0.01
    push!(koz_min, [middle[1]-a,  middle[2]+b,     -c])
    push!(koz_max, [middle[1]+a,  worldAABBmax[2],  c])
    push!(koz_min, [middle[1]-a,  worldAABBmin[2], -c])
    push!(koz_max, [middle[1]+a,  middle[2]-b,      c])

    # --------------------------
  end


  for i in 1:length(koz_min)
    push!(keepout_zones, HyperRectangle(koz_min[i]..., (koz_max[i]-koz_min[i])...))
  end
  obstacle_set = Vector{GeometryTypes.GeometryPrimitive}(0)


  # -----------------------------------
  # CoRL- Additional start / end  zones

  # Safe region
  safe_region_type = "rectangle"
  nb_edges = 4
  rght_safe_x = 1.35
  left_safe_x = 0.3
  top_safe_y  = 2.4
  btm_safe_y  = 0.3
  safe_region_coeffs = fill(Float64[], nb_edges)
  safe_region_coeffs[1] = [1.,  0., 0., -rght_safe_x]
  safe_region_coeffs[2] = [-1., 0., 0.,  left_safe_x]
  safe_region_coeffs[3] = [0.,  1., 0., -top_safe_y ]
  safe_region_coeffs[4] = [0., -1., 0.,  btm_safe_y ]
  safe_region_smallest_width = min(rght_safe_x-left_safe_x, top_safe_y-btm_safe_y)
  safe_region_height = top_safe_y  - btm_safe_y
  safe_region_width  = rght_safe_x - left_safe_x

  #  Goal / end region
  # goal_region_type = "ball"
  # goal_region_position = [0.85; 0.8; 0.]
  # goal_region_radius   = 0.6
  # -----------------------------------


  return Table{T}(worldAABBmin, worldAABBmax, keepin_zones, keepout_zones, obstacle_set,
                  safe_region_type, safe_region_coeffs, 
                  safe_region_smallest_width, safe_region_height, safe_region_width)
end
