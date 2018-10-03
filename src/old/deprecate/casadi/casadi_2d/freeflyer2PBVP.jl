using Convex, Gurobi
using PyPlot, HDF5

include("../../robot.jl")
include("../../sampling.jl")

function generate_samples()
  rb = astrobee()
  Xmin1,Xmin2,Xmax = -0.5,-0.75,0.75
  Ymin,Ymax = -0.75,0.75
  rb = astrobee()

  n_samples = 500

  # generate samples from ([-0.75,0.75], [-0.75,0.75])
  Vfull = generateHaltonSamples(4,2*n_samples,[-(Xmax-Xmin2)/2,-(Ymax-Ymin)/2,-rb.hard_limit_vel/sqrt(2),-rb.hard_limit_vel/sqrt(2)])

  # prune samples ([-0.5,0.75], [-0.75,0.75])
  valid_idx = find(Vfull[1,:] .> Xmin1)
  V = Vfull[:,valid_idx[1:n_samples]]

  # manually add desired points
  V[:,1] = zeros(4)
  n_grid = 10
  grid = collect(Iterators.product(linspace(Xmin1,Xmax,n_grid), linspace(Ymin,Ymax,n_grid)))
  X = collect(Iterators.flatten(map(x->x[1],grid))) 
  Y = collect(Iterators.flatten(map(x->x[2],grid))) 
  # V[1,2:1+n_grid^2] = X
  # V[2,2:1+n_grid^2] = Y 

  h5open("V.h5", "w") do file
    write(file, "V", V)
  end
end

function setup()
  n,m = 4,2
  tf, N = 25., 51
  dt = tf/(N-1)

  rb = astrobee()

  X = Variable(n,N)
  U = Variable(m,N-1)

  constraints = Convex.Constraint[]

  Ak = kron([1 dt; 0 1], eye(Int(n/2)))
  Bk = kron([0.5*dt^2; dt], eye(Int(n/2)))

  for k in 1:N
    constraints += norm(X[3:4,k],2) <= rb.hard_limit_vel 
  end

  for k in 1:N-1
    constraints += X[:,k+1] == Ak*X[:,k]+Bk*U[:,k]
    constraints += norm(U[:,k],2) <= rb.hard_limit_accel
  end

  constraints += X[:,1] == zeros(n) 
  constraints += X[:,end] == zeros(n)

  return constraints, X, U
end

function update{T}(Xi::Vector{T},Xf::Vector{T},constraints)
  constraints[end-1] = X[:,1] == Xi
  constraints[end] = X[:,end] == Xf

  Jm = quadform(U,eye(m))
  prob = minimize(Jm,constraints)
  @time Convex.solve!(prob,GurobiSolver())
end

function main()
  set_default_solver(GurobiSolver())
  env = Gurobi.Env()
  setparam!(env,"OutputFlag",0)

  constraints,X,U = setup()
  n,m = size(X,1), size(U,1)

  V = h5read("V.h5", "V")

  n_samples = size(V,2)
  max_dist = 0.35

  fn = "costs_$(Int(100*max_dist))cm.h5"

  costs = Inf*ones(n_samples^2)

  t_cum = 0.
  @time for i = 327:n_samples
    tic()
    for j  = 1:n_samples
      dist = norm(V[1:2,i]-V[1:2,j])

      if i==j
        continue
        println("Not calculating 2PBVP for ($i,$j)")
      elseif dist > max_dist
        println("Not calculating 2PBVP for ($i,$j)")
        continue
      end

      constraints[end-1] = X[:,1] == V[:,i] 
      constraints[end] = X[:,end] == V[:,j]

      Jm = quadform(U,eye(m))
      prob = minimize(Jm,constraints)
      @time Convex.solve!(prob)

      if prob.status == :Optimal
        costs[n_samples*(i-1)+j] = prob.optval
        println("Optimal soln found for ($i,$j)")
      end
      println("Optimal soln not found for ($i,$j)")

      gc() 
    end
    
    h5open(fn, "w") do file
      write(file, "costs", costs)
    end

    t_tot = toc()
    t_cum += t_tot
    println("Finished calculating TPBVP's for node $i")
    println("TPBVPs time for node $i was $t_tot")
    println("Time so far $t_cum")
  end
  h5open(fn, "w") do file
    write(file, "costs", costs)
  end
  print("Done")
end
