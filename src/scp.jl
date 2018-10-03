
include("scp/scp_gusto.jl")
include("scp/scp_trajopt.jl")
include("scp/scp_mao.jl")

function convergence_metric(traj::Trajectory, traj_prev::Trajectory, SCPP::SCPProblem)
  # normalized maximum relative error between iterations
  max_num, max_den = -Inf, -Inf
  for k in 1:SCPP.N
    val = norm(traj.X[:,k]-traj_prev.X[:,k])
    max_num = val > max_num ? val : max_num

    val = norm(traj.X[:,k])
    max_den = val > max_den ? val : max_den
  end
  return max_num/max_den
end
