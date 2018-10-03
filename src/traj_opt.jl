
export solve_SCPshooting!, solve_SCP!

function solve_SCPshooting!(TOS::TrajectoryOptimizationSolution, TOP::TrajectoryOptimizationProblem, solve_method!, init_method, solver="Mosek"; kwarg...)
	robot = TOP.PD.robot
	model = TOP.PD.model

	SCPP = SCPProblem(TOP)

	# Create initialization trajectory
	traj_init = init_method(TOP)

	# Run SCP with initialization trajectory
	SCPS = SCPSolution(SCPP, traj_init)
	SP = ShootingProblem(TOP, SCPS)
	SS = ShootingSolution()

	solve_method!(SCPS, SCPP, solver, 1; kwarg...)

	# Until shooting method succeeds or maximum SCP iterations is reached
	while (!SCPS.converged)
		# Attempt shooting method
		SP = ShootingProblem(TOP, SCPS)
		solve!(SS, SP)
		
		# If shooting converged, exit
		SS.converged ? break : nothing

		# Run another iteration of SCP
		solve_method!(SCPS, SCPP, solver, 1; kwarg...)
	end

	# If a solution has been found and there is no time remaining for refinement, return solution
	TOS.traj, TOS.SCPS, TOS.SS = SCPS.traj, SCPS, SS
end

function solve_SCP!(TOS::TrajectoryOptimizationSolution, TOP::TrajectoryOptimizationProblem, solve_method!, init_method, solver="Mosek"; kwarg...)
	robot = TOP.PD.robot
	model = TOP.PD.model

	SCPP = SCPProblem(TOP)

	# Create initialization trajectory
	traj_init = init_method(TOP)

	# Run SCP with initialization trajectory
	SCPS = SCPSolution(SCPP, traj_init)
	TOS.traj, TOS.SCPS = SCPS.traj, SCPS
	solve_method!(SCPS, SCPP, solver; kwarg...)
end
