using GeometryTypes

include("mpproblem.jl")

function is_free{T}(x::Vector{T},zones::Vector{HyperRectangle})
	# https://math.stackexchange.com/questions/1472049/check-if-a-point-is-inside-a-rectangular-shaped-area-3d
  for zone in zones
    lo,hi = zone.origin, zone.origin+zone.widths
    v = x - lo
    in_collision = true 

    for i in 1:3
      u = zeros(T,3)
      u[i] = hi[i]-lo[i]
      in_collision &= (0 <= dot(v,u) <= dot(u,u))
    end
    in_collision && return false 
  end
  return true
end

function is_free_edge{T}(Xedge::Matrix{T},P::MPProblem{T})
  check_idx = [Int(ceil(i*size(Xedge,2)/P.n_waypoints+1)) for i in 0:P.n_waypoints-1]
  for (k,idx) in enumerate(check_idx)
    if !is_free(Xedge[1:3,idx],P.world.keepout_zones) || 
        is_free(Xedge[1:3,idx],P.world.keepin_zones)
      return false
    end
  end
  return true
end
