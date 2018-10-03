include("quat_functions.jl")

function tpbvp_cost{T<:AbstractFloat}(tau::T,xi::Vector{T},xf::Vector{T},R::Matrix{T})
  return tau + ((6*R[1,1]*((xf[4]) - (xi[4])))/tau^2 + (12*R[1,1]*((xi[1]) - (xf[1]) + (tau)*(xi[4])))/tau^3)*(xi[1] - xf[1] + tau*xi[4]) + ((6*R[2,2]*((xf[5]) - (xi[5])))/tau^2 + (12*R[2,2]*((xi[2]) - (xf[2]) + (tau)*(xi[5])))/tau^3)*(xi[2] - xf[2] + tau*xi[5]) + ((6*R[3,3]*((xf[6]) - (xi[6])))/tau^2 + (12*R[3,3]*((xi[3]) - (xf[3]) + (tau)*(xi[6])))/tau^3)*(xi[3] - xf[3] + tau*xi[6]) + ((4*R[1,1]*((xf[4]) - (xi[4])))/tau + (6*R[1,1]*((xi[1]) - (xf[1]) + (tau)*(xi[4])))/tau^2)*(xf[4] - xi[4]) + ((4*R[2,2]*((xf[5]) - (xi[5])))/tau + (6*R[2,2]*((xi[2]) - (xf[2]) + (tau)*(xi[5])))/tau^2)*(xf[5] - xi[5]) + ((4*R[3,3]*((xf[6]) - (xi[6])))/tau + (6*R[3,3]*((xi[3]) - (xf[3]) + (tau)*(xi[6])))/tau^2)*(xf[6] - xi[6])
end


function dcost{T<:AbstractFloat}(tau::T,xi::Vector{T},xf::Vector{T},R::Matrix{T})
  dtau = 0.000001
	return (tpbvp_cost(tau+dtau/2,xi,xf,R) - tpbvp_cost(tau-dtau/2,xi,xf,R))/dtau
end


function toptBisection{T}(x0::Vector{T},x1::Vector{T},tmin::T,tmax::T,R::Matrix{T})
  TOL = 0.0000001
	tu = tmax
  
  if (dcost(tmax,x0,x1,R) < 0) 
    return tmax
  end 

  tl = 0.01
  tl = tmin
  while (dcost(tl,x0,x1,R) > 0) 
    tl = tl/2
  end

  topt = 0
  dcval = 1
  maxItrs = 80
  itr = 0

  while (abs(dcval) > TOL && itr < maxItrs) 
    topt = (tu+tl)/2
    dcval = dcost(topt,x0,x1,R)
    if (dcval > 0)
      tu = topt
    else 
      tl = topt
    end

    itr += 1
  end
  return topt
end


function pathPoint{T<:AbstractFloat}(t::T,tau::T,xi::Vector{T},xf::Vector{T},R::Matrix{T})
  # t: time to query
  # tau: final time
  x = similar(xi) 
  x[1] = xf[1] + xf[4]*(t - tau) + (((4*R[1,1]*(xf[4] - xi[4]))/tau + (6*R[1,1]*(xi[1] - xf[1] + tau*xi[4]))/tau^2)*(t - tau)^2)/(2*R[1,1]) + (((6*R[1,1]*(xf[4] - xi[4]))/tau^2 + (12*R[1,1]*(xi[1] - xf[1] + tau*xi[4]))/tau^3)*(t - tau)^3)/(6*R[1,1])
  x[2] = xf[2] + xf[5]*(t - tau) + (((4*R[2,2]*(xf[5] - xi[5]))/tau + (6*R[2,2]*(xi[2] - xf[2] + tau*xi[5]))/tau^2)*(t - tau)^2)/(2*R[2,2]) + (((6*R[2,2]*(xf[5] - xi[5]))/tau^2 + (12*R[2,2]*(xi[2] - xf[2] + tau*xi[5]))/tau^3)*(t - tau)^3)/(6*R[2,2])
  x[3] = xf[3] + xf[6]*(t - tau) + (((4*R[3,3]*(xf[6] - xi[6]))/tau + (6*R[3,3]*(xi[3] - xf[3] + tau*xi[6]))/tau^2)*(t - tau)^2)/(2*R[3,3]) + (((6*R[3,3]*(xf[6] - xi[6]))/tau^2 + (12*R[3,3]*(xi[3] - xf[3] + tau*xi[6]))/tau^3)*(t - tau)^3)/(6*R[3,3])
  x[4] = xf[4] + (((4*R[1,1]*(xf[4] - xi[4]))/tau + (6*R[1,1]*(xi[1] - xf[1] + tau*xi[4]))/tau^2)*(t - tau))/R[1,1] + (((6*R[1,1]*(xf[4] - xi[4]))/tau^2 + (12*R[1,1]*(xi[1] - xf[1] + tau*xi[4]))/tau^3)*(t - tau)^2)/(2*R[1,1])
  x[5] = xf[5] + (((4*R[2,2]*(xf[5] - xi[5]))/tau + (6*R[2,2]*(xi[2] - xf[2] + tau*xi[5]))/tau^2)*(t - tau))/R[2,2] + (((6*R[2,2]*(xf[5] - xi[5]))/tau^2 + (12*R[2,2]*(xi[2] - xf[2] + tau*xi[5]))/tau^3)*(t - tau)^2)/(2*R[2,2])
  x[6] = xf[6] + (((4*R[3,3]*(xf[6] - xi[6]))/tau + (6*R[3,3]*(xi[3] - xf[3] + tau*xi[6]))/tau^2)*(t - tau))/R[3,3] + (((6*R[3,3]*(xf[6] - xi[6]))/tau^2 + (12*R[3,3]*(xi[3] - xf[3] + tau*xi[6]))/tau^3)*(t - tau)^2)/(2*R[3,3])

  # TODO(acauligi): Check to make sure slerp function is correct and matches Eigen library/Ames convention
  # x[7:10] = quat_interp(xi[7:10],xf[7:10],t/tau)
  return x
end
