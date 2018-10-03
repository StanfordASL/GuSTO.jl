using HDF5
include("dynamics.jl")

function t_to_tau{T}(t::Vector{T},tau0::T,tauf::T)
  """ Converts from the interval [-1,1] to the interval [tau0,tauf] """
  N = length(t)-1
  tau = zeros(T,N+1)
  for i in 1:N+1 
    tau[i] = 0.5*((tauf-tau0)*t[i] + (tauf + tau0))
  end
  return tau
end

function tau_to_t{T}(tau::Vector{T})
  """ Converts from the interval [tau0,tauf] to the interval [-1,1] """
  N = length(tau)-1
  t = zeros(T,N+1)
  tau0 = tau[0]
  tauf = tau[end]
  for i in 1:N+1
    t[i] = (2.0*tau[i] - (tauf-tau0))/(tauf-tau0)
  end
  return t
end

function chebyshev_nodes(N::Int)
  """ Returns N+1 CGL points on domain [1,-1] for the extrema
  of the Nth order Chebyshev polynomial. """
  x = zeros(N+1)
  for k in 0:N
    x[k+1] = cos(k*pi/N)
  end
  return x[end:-1:1]
end

function polyApproximation{T}(xvals::Vector{T},yvals::Vector{T})
  function LagrangeInterpolant{T}(x::T)
    N = length(xvals)
    LagrangePolynomials = ones(T,N)
    for i in 1:N  
      for j in [1:i-1;i+1:N]
        LagrangePolynomials[i] = LagrangePolynomials[i].*(x-xvals[j])./(xvals[i]-xvals[j])
      end
    end
    output = sum(LagrangePolynomials.*yvals)
  end
  return LagrangeInterpolant
end

# function lagrange{T}(ts::Vector{T},slice::Vector{T})
#   N = length(slice)-1
#   P = 1.
#   for j in 1:N+1
#     P_here = slice[j] 
#     for k in 1:N+1
#       if j!=k 
#         P_here *= 1/(ts[j]-ts[k])*Poly([-ts[k],1.], :r)
#       end
#     end
#     P += P_here
#   end
#   return P
# end
