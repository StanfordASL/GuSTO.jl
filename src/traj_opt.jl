
export solve_SCPshooting!, solve_SCP!

function solve_SCPshooting!(TOS::TrajectoryOptimizationSolution, TOP::TrajectoryOptimizationProblem, solve_method!, init_method, solver="Mosek"; kwarg...)
	robot = TOP.PD.robot
	model = TOP.PD.model

	SCPP = SCPProblem(TOP)

	# Create initialization trajectory
	traj_init = init_method(TOP)

	# Run SCP with initialization trajectory
	TOS.SCPS = SCPSolution(SCPP, traj_init)
	SCPS = TOS.SCPS
	SP = ShootingProblem(TOP, SCPS)
	TOS.SS = ShootingSolution(SP, deepcopy(traj_init))
	SS = TOS.SS

	solve_method!(SCPS, SCPP, solver, 1; kwarg...)
	push!(SS.J_true, SCPS.J_true[1])

	# Until shooting method succeeds or maximum SCP iterations is reached
	ss_sol = nothing
	while (!SCPS.converged)
		# Attempt shooting method
		SP = ShootingProblem(TOP, SCPS)
		ss_sol = solve!(SS, SP)
		
		# If successful shooting runs have converged, exit
		conv_iter_spread = 3
		if SCPS.iterations > conv_iter_spread && sum(SS.convergence_measure[end-conv_iter_spread+1:end]) <= SCPS.param.convergence_threshold
			SS.converged = true
			# TODO: Check inequality constraints
			copy!(TOS.traj, SS.traj)
			TOS.total_time = SCPS.total_time + sum(SS.iter_elapsed_times)
			return
		end
		# Run another iteration of SCP
		solve_method!(SCPS, SCPP, solver, 1; kwarg...)
	end

	copy!(TOS.traj, SCPS.traj)
	TOS.total_time = SCPS.total_time + sum(SS.iter_elapsed_times)
end

function solve_SCP!(TOS::TrajectoryOptimizationSolution, TOP::TrajectoryOptimizationProblem, solve_method!, init_method, solver="Mosek"; max_iter=50, force=false, kwarg...)
	robot = TOP.PD.robot
	model = TOP.PD.model

	SCPP = SCPProblem(TOP)

	# Create initialization trajectory
	traj_init = init_method(TOP)

	# Run SCP with initialization trajectory
	SCPS = SCPSolution(SCPP, traj_init)
	TOS.traj, TOS.SCPS = SCPS.traj, SCPS
	solve_method!(SCPS, SCPP, solver, max_iter, force; kwarg...)
end

function solve_SCP!(TOS::TrajectoryOptimizationSolution, TOP::TrajectoryOptimizationProblem, solve_method!, traj_init::Trajectory, solver="Mosek"; max_iter=50, force=false, kwarg...)
	robot = TOP.PD.robot
	model = TOP.PD.model

	SCPP = SCPProblem(TOP)
	
	# Run SCP with initialization trajectory
	SCPS = SCPSolution(SCPP, traj_init)
	TOS.traj, TOS.SCPS = SCPS.traj, SCPS
	solve_method!(SCPS, SCPP, solver, max_iter, force; kwarg...)
end

function convergence_metric(traj::Trajectory, traj_prev::Trajectory, OAP::OptAlgorithmProblem)
  # normalized maximum relative error between iterations
  max_num, max_den = -Inf, -Inf
  for k in 1:OAP.N
    val = norm(traj.X[:,k]-traj_prev.X[:,k])
    max_num = val > max_num ? val : max_num

    val = norm(traj.X[:,k])
    max_den = val > max_den ? val : max_den
  end
  return max_num/max_den
end