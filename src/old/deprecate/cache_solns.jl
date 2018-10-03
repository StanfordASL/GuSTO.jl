using HDF5

include("load_stuff.jl")

# function writeSolution{T<:AbstractString}(status_fn::T,edges_fn::T,tpbvps_fn::T,times::Vector,tran_costs::Vector,rot_costs::Vector,tpbvps::Matrix, P::MPProblem)
#   # Edge data
#   x_dim           = h5read(edges_fn, "x_dim")
#   n_samples       = h5read(edges_fn, "n_samples")
#   n_waypoints     = h5read(edges_fn, "n_waypoints")
#   R               = h5read(edges_fn, "R")
#   iteration       = h5read(edges_fn, "iteration")
#   old_times       = h5read(edges_fn, "times")
#   old_tran_costs  = h5read(edges_fn, "tran_costs")
#   old_rot_costs   = h5read(edges_fn, "rot_costs")
# 
#   times       = convert(Array{Float64}, times)
#   tran_costs  = convert(Array{Float64}, tran_costs)
#   rot_costs   = convert(Array{Float64}, rot_costs)
# 
#   h5open(edges_fn, "w") do file
#     write(file, "x_dim", x_dim)
#     write(file, "n_samples", n_samples)
#     write(file, "n_waypoints", n_waypoints)
#     write(file, "iteration", iteration+1)
#     write(file, "times", [old_times; times])
#     write(file, "R", R)
#     write(file, "tran_costs", [old_tran_costs; tran_costs])
#     write(file, "rot_costs", [old_rot_costs; rot_costs])
#   end
# 
#   old_times = []
#   old_rot_costs = []
#   old_tran_costs = []
#   gc()
#   
#   # # TPBVP data
#   # (x,y,z) = size(tpbvps)
#   # new_tpbvps = reshape(tpbvps, (x*z,y))
#   # old_tpbvps = h5read(tpbvps_fn, "tpbvps")
#   # h5open(tpbvps_fn, "w") do file
#   #   write(file, "tpbvps", [old_tpbvps; new_tpbvps])
#   # end
# 
#   old_tpbvps = []
#   new_tpbvps = []
#   gc()
# end

use3D = true
P = MPProblem()

# if use3D
#   env = P.world
#   rb = P.robot
#   sample_max_speed = -rb.hard_limit_vel/sqrt(3) # want norm(v) <= max_speed
#   dim = 6     # sample over R^6
# 
#   rv_samples = generateHaltonSamples(dim, P.n_samples, [-(env.dim_x-rb.r); -(env.dim_y-rb.r); -(env.dim_z-rb.r); sample_max_speed; sample_max_speed; sample_max_speed])
#   rv_samples[1:3,:] -= 1
# 
#   # Rotational planning samples
#   q0 = [0. 0. 0. 1.]'
#   q_samples = [q0 yershova_quaternions()[:, 1:sp.n_samples-1]]
# 
#   samples = [valid_rv_samples; q_samples]
# else
#   env = table()
#   # Generate samples
#   sample_max_speed = -rb.hard_limit_vel/sqrt(2) # want norm(v) <= max_speed
#   dim = 7     # sample over SE(2)
# 
#   se2_samples = generateHaltonSamples(dim, 10*sp.n_samples, [-(0.5*env.dim_x-rb.r); -(0.5*env.dim_y-rb.r); 0; 
#                                                           -sample_max_speed; -sample_max_speed; 0; 2*pi])
#   samples = zeros(10,sp.n_samples)
# 
#   for idx in 1:sp.n_samples
#     samples[1:3,idx] = se2_samples[1:3,idx]
#     samples[4:6,idx] = se2_samples[4:6,idx]
#     samples[7:10,idx] = vec2quat(env.up,se2_samples[7,idx])
#   end
# 
#   # Initial node
#   samples[1:2,1] = -0.5*ones(2)
#   samples[4:6,1] = zeros(3)
#   samples[7:10,1] = vec2quat(tb.up,0.)
# end

stringtime = string(string(Dates.today()), "_", string(Dates.hour(now()), "_", string(Dates.minute(now()))))
tpbvps_fn = string("tpbvps_", stringtime, ".h5")
# samples_fn = string("samples_set_", stringtime, ".h5")
# edges_fn = string("edges_", stringtime, ".h5")
# status_fn = string("status_", stringtime, ".h5")

# Seed files
# h5open(samples_fn, "w") do file
#   write(file, "samples", P.V)
# end

# h5open(edges_fn, "w") do file
#   write(file, "x_dim", P.x_dim)
#   write(file, "n_samples", P.n_samples)
#   write(file, "n_waypoints", P.n_waypoints)
#   write(file, "iteration", 0)
#   write(file, "times", Vector{Float64}())
#   write(file, "R", P.R)
#   write(file, "tran_costs",Vector{Float64}())
#   write(file, "rot_costs", Vector{Float64}())
# end

# h5open(tpbvps_fn, "w") do file
#   write(file, "tpbvps", Array{Float64}(1,P.n_waypoints))
# end

# Generating TPBVP's
times = Float64(Inf)*ones(Float64,P.n_samples^2) 
tran_costs = Float64(Inf)*ones(Float64,P.n_samples^2) 
rot_costs = Float64(Inf)*ones(Float64,P.n_samples^2) 

idx = 1
@time for i in 1:P.n_samples
  tic()
  X0 = P.V[:,i]
  # Xtpbvp = zeros(P.x_dim,P.n_waypoints,P.n_samples)
  for j in 1:P.n_samples
    if i != j
      Xf = P.V[:,j]
      topt = toptBisection(X0,Xf,P.t_range[1],P.t_range[2],P.R)
      J = tpbvp_cost(topt,X0,Xf,P.R)
      times[idx] = topt 
      tran_costs[idx] = J

      dot = X0[7:10]'*Xf[7:10]
      if (abs(dot) > 1.)
        # check against numerical error
        dot = sign(dot)
      end
      rot_costs[idx] = acos(dot)
      # rot_costs[idx] = quat2angle(quat_error(X0[7:10], Xf[7:10]))
      # rot_costs[idx] = norm(quat_log(quat_error(X0[7:10],Xf[7:10])))
      # [Xtpbvp[:,m,j] = pathPoint(t,topt,X0,Xf,P.R) for (m,t) in enumerate(linspace(0,topt,P.n_waypoints))]
    end
    idx+=1
  end

  # writeSolution(status_fn, edges_fn, tpbvps_fn, times, tran_costs, rot_costs, Xtpbvp, P)
  t_total = toc()
  println("Finished calculating TPBVP's for node $i")
  println("TPBVPs time for node $i was $t_total")
end

h5open(tpbvps_fn, "w") do file
  write(file, "V", P.V)
  write(file, "x_dim", P.x_dim)
  write(file, "n_samples", P.n_samples)
  write(file, "n_waypoints", P.n_waypoints)
  write(file, "R", P.R)
  write(file, "times", times)
  write(file, "tran_costs",tran_costs)
  write(file, "rot_costs", rot_costs)
end
