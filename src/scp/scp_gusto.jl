
export solve_gusto_jump!#, solve_gusto_jump_nonlinear!

mutable struct SCPParam_GuSTO <: SCPParamSpecial
	Δ0 								# trust region size 
	ω0         				# exact penalty violation
	ω_max
	ε
	ρ0 
	ρ1
  β_succ		   				# trust-region growth factors
  β_fail 	
  γ_fail

  Δ_vec::Vector
  ω_vec::Vector
  ρ_vec::Vector           # trust region ratios
  trust_region_satisfied_vec::Vector  # Bool flags
  convex_ineq_satisfied_vec::Vector  # Bool flags
end

function SCPParam_GuSTO(Δ0, ω0, ω_max, ε, ρ0, ρ1, β_succ, β_fail, γ_fail)
	SCPParam_GuSTO(Δ0, ω0, ω_max, ε, ρ0, ρ1, β_succ, β_fail, γ_fail, [Δ0], [ω0], [0.], [false], [false])
end

##########
# General
##########

# Define this for a dynamics model if model parameters need to be updated
# at each SCP iteration
update_model_params!(SCPP::SCPProblem, traj_prev::Trajectory) = nothing

function trust_region_satisfied_gusto(traj::Trajectory, traj_prev::Trajectory, SCPP::SCPProblem)
  # checks for satisfaction of state trust region constraint
  Δ = SCPP.param.alg.Δ_vec[end]
  max_val = -Inf

  for k in 1:SCPP.N
    val = norm(traj.X[:,k]-traj_prev.X[:,k])^2
    max_val = val > max_val ? val : max_val
  end
  return max_val-Δ <= 0
end

########
# JuMP
########
function solve_gusto_jump!(SCPS::SCPSolution, SCPP::SCPProblem, solver="Ipopt", max_iter=30, force=false; kwarg...)
	# Solves a sequential convex programming problem
	# Inputs:
	#		SCPS - SCP solution structure, including an initial trajectory
	#		SCPP - SCP problem definition structure
	#   max_iter - Set to enforce a maximum number of iterations (0 = unlimited)
	#   force - Set true to force further refinement iterations despite convergence	

	N = SCPP.N
	param, model = SCPP.param, SCPP.PD.model

	!isdefined(param, :alg) ? param.alg = SCPParam_GuSTO(model) : nothing
	Δ0, ω0, ω_max, ρ0, ρ1 = param.alg.Δ0, param.alg.ω0, param.alg.ω_max, param.alg.ρ0, param.alg.ρ1
	β_succ, β_fail, γ_fail = param.alg.β_succ, param.alg.β_fail, param.alg.γ_fail
	Δ_vec, ω_vec, ρ_vec = param.alg.Δ_vec, param.alg.ω_vec, param.alg.ρ_vec
	trust_region_satisfied_vec, convex_ineq_satisfied_vec = param.alg.trust_region_satisfied_vec, param.alg.convex_ineq_satisfied_vec

	# TODO: Modify constraints rather than creating new problem
	iter_cap = SCPS.iterations + max_iter
	SCPV = SCPVariables{JuMP.VariableRef, Array{JuMP.VariableRef}}()
	SCPS.SCPC = SCPConstraints(SCPP)
	SCPC = SCPS.SCPC

	initialize_model_params!(SCPP, SCPS.traj)
	push!(SCPS.J_true, cost_true(SCPS.traj, SCPS.traj, SCPP))
	push!(SCPS.J_full, SCPS.J_true[end])
	push!(ρ_vec, trust_region_ratio_gusto(SCPS.traj, SCPS.traj, SCPP))
	param.obstacle_toggle_distance = Δ_vec[end]/8 + model.clearance # TODO: Generalize clearance

	while (SCPS.iterations < iter_cap)
		time_start = time_ns()

		# TODO: Avoid reconstructing the problem from scratch
		if solver == "Mosek"
			SCPS.solver_model = Model(with_optimizer(Mosek.Optimizer; kwarg...))
		elseif solver == "Gurobi"
			SCPS.solver_model = Model(with_optimizer(Gurobi.Optimizer; kwarg...))
		elseif solver == "Ipopt"
			SCPS.solver_model = Model(with_optimizer(Ipopt.Optimizer; kwarg...))
		elseif solver == "SCS"
			SCPS.solver_model = Model(with_optimizer(SCS.Optimizer; kwarg...))
		else 	# Default solver
			SCPS.solver_model = Model(with_optimizer(Ipopt.Optimizer; kwarg...))
		end

		# Set up, solve problem
		update_model_params!(SCPP, SCPS.traj)
		add_variables_jump!(SCPS, SCPV, SCPP)
		add_constraints_gusto_jump!(SCPS, SCPV, SCPC, SCPP)
		add_objective_gusto_jump!(SCPS, SCPV, SCPC, SCPP)

		set_start_value.(SCPV.X, SCPS.traj.X)
		set_start_value.(SCPV.U, SCPS.traj.U)
		set_start_value(SCPV.Tf, SCPP.tf_guess)

		JuMP.optimize!(SCPS.solver_model)

		push!(SCPS.solver_status, JuMP.termination_status(SCPS.solver_model))
		if SCPS.solver_status[end] ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED)
		  # warn("GuSTO SCP iteration failed to find an optimal solution")
		  push!(SCPS.iter_elapsed_times, (time_ns() - time_start)/10^9) 
		  return
		end

		# Recover solution
		new_traj = Trajectory(JuMP.value.(SCPV.X), JuMP.value.(SCPV.U), JuMP.value.(SCPV.Tf))
		push!(SCPS.convergence_measure, convergence_metric(new_traj, SCPS.traj, SCPP))
		push!(SCPS.J_full, JuMP.objective_value(SCPS.solver_model))
		SCPS.dual = get_dual_jump(SCPC, SCPP)
		
		# Check trust regions
		push!(trust_region_satisfied_vec, trust_region_satisfied_gusto(new_traj, SCPS.traj, SCPP))
		push!(convex_ineq_satisfied_vec, convex_ineq_satisfied_gusto_jump(new_traj, new_traj, SCPC, SCPP))

		if trust_region_satisfied_vec[end]
			push!(ρ_vec, trust_region_ratio_gusto(new_traj, SCPS.traj, SCPP))			
			if ρ_vec[end] > ρ1
				push!(SCPS.scp_status, :InaccurateModel)
				push!(SCPS.accept_solution, false)
				push!(Δ_vec, β_fail*Δ_vec[end])
				push!(ω_vec, ω_vec[end])
			else
				push!(SCPS.accept_solution, true)
				ρ_vec[end] < ρ0 ? push!(Δ_vec, min(β_succ*Δ_vec[end], Δ0)) : push!(Δ_vec, Δ_vec[end])
				if !convex_ineq_satisfied_vec[end]
					push!(SCPS.scp_status, :ViolatesConstraints)
					push!(ω_vec, γ_fail*ω_vec[end])
				else
					push!(SCPS.scp_status, :OK)
					# push!(ω_vec, ω0)
					push!(ω_vec, ω_vec[end])
				end
			end
		else
			push!(SCPS.scp_status, :TrustRegionViolated)
			push!(SCPS.accept_solution, false)
			push!(Δ_vec, Δ_vec[end])
			push!(ω_vec, γ_fail*ω_vec[end])
		end

		if SCPS.accept_solution[end]
			push!(SCPS.J_true, cost_true(new_traj, SCPS.traj, SCPP))
			copy!(SCPS.traj, new_traj)	# TODO(ambyld): Maybe deepcopy not needed?
		else
			push!(SCPS.J_true, SCPS.J_true[end])
		end

		param.obstacle_toggle_distance = Δ_vec[end]/8 + model.clearance

		iter_elapsed_time = (time_ns() - time_start)/10^9
		push!(SCPS.iter_elapsed_times, iter_elapsed_time)
		SCPS.total_time += iter_elapsed_time
		SCPS.iterations += 1

		if ω_vec[end] > ω_max
			warn("GuSTO SCP omegamax exceeded")
			break
		end
		!SCPS.accept_solution[end] ? continue : nothing

		conv_iter_spread = 2
		if SCPS.iterations > conv_iter_spread && sum(SCPS.convergence_measure[end-conv_iter_spread+1:end]) <= param.convergence_threshold
			SCPS.converged = true
			convex_ineq_satisfied_vec[end] && (SCPS.successful = true)
			force ? continue : break
		end
	end
end

function add_variables_jump!(SCPS::SCPSolution, SCPV::SCPVariables, SCPP::SCPProblem)
	solver_model = SCPS.solver_model
	model = SCPP.PD.model
  x_dim, u_dim, N = model.x_dim, model.u_dim, SCPP.N

	@variable(solver_model, X[1:x_dim, 1:N])
	@variable(solver_model, U[1:u_dim, 1:N])
	@variable(solver_model, Tf)

	SCPP.param.fixed_final_time ? JuMP.fix(Tf, SCPP.tf_guess) : nothing

	SCPV.X, SCPV.U, SCPV.Tf = X, U, Tf
end

function add_constraints_gusto_jump!(SCPS::SCPSolution, SCPV::SCPVariables, SCPC::SCPConstraints, SCPP::SCPProblem)
	solver_model, traj_prev = SCPS.solver_model, SCPS.traj

	update_model_params!(SCPP, traj_prev)

	for cc_list in values(SCPC.dynamics)
		for cc in cc_list
			if cc.dimtype == :scalar
				cc.con_reference = @constraint(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)], cc.func(SCPV, traj_prev, SCPP, k, i...) == 0)
			elseif cc.dimtype == :array
				cc.con_reference = @constraint(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)], cc.func(SCPV, traj_prev, SCPP, k, i...) .== 0)
			end
		end
	end

	for cc_list in values(SCPC.state_init_eq)
		for cc in cc_list
			if cc.dimtype == :scalar
				cc.con_reference = @constraint(solver_model, [i=Iterators.product(cc.ind_other...)], cc.func(SCPV, traj_prev, SCPP, 1, i...) == 0)
			elseif cc.dimtype == :array
				cc.con_reference = @constraint(solver_model, [i=Iterators.product(cc.ind_other...)], cc.func(SCPV, traj_prev, SCPP, 1, i...) .== 0)
			end
		end
	end
	
  for cc_list in values(SCPC.convex_control_ineq)
  	for cc in cc_list
			if cc.dimtype == :scalar
				cc.con_reference = @constraint(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)], cc.func(SCPV, traj_prev, SCPP, k, i...) <= 0)
			elseif cc.dimtype == :array
				cc.con_reference = @constraint(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)], cc.func(SCPV, traj_prev, SCPP, k, i...) .<= 0)
			end
		end
	end

	for cc_list in values(merge(SCPC.convex_state_boundary_condition_eq, SCPC.nonconvex_state_boundary_condition_convexified_eq))
		for cc in cc_list
			if cc.dimtype == :scalar
				cc.con_reference = @constraint(solver_model, [k=cc.params.k_timestep,i=Iterators.product(cc.ind_other...)], cc.func(SCPV, traj_prev, SCPP, cc.params, i...) == 0)
			elseif cc.dimtype == :array
				cc.con_reference = @constraint(solver_model, [k=cc.params.k_timestep,i=Iterators.product(cc.ind_other...)], cc.func(SCPV, traj_prev, SCPP, cc.params, i...) .== 0)
			end
		end
	end

	for cc_list in values(merge(SCPC.convex_state_boundary_condition_ineq, SCPC.nonconvex_state_boundary_condition_convexified_ineq))
		for cc in cc_list
			if cc.dimtype == :scalar
				cc.con_reference = @constraint(solver_model, [k=cc.params.k_timestep,i=Iterators.product(cc.ind_other...)], cc.func(SCPV, traj_prev, SCPP, cc.params, i...) <= 0)
			elseif cc.dimtype == :array
				cc.con_reference = @constraint(solver_model, [k=cc.params.k_timestep,i=Iterators.product(cc.ind_other...)], cc.func(SCPV, traj_prev, SCPP, cc.params, i...) .<= 0)
			end
		end
	end

	# TODO: Store constraint reference to this
	if !SCPP.param.fixed_final_time
		@constraint(solver_model, SCPV.Tf >= 0.1)
	end
end

function add_objective_gusto_jump!(SCPS::SCPSolution, SCPV::SCPVariables, SCPC::SCPConstraints, SCPP::SCPProblem)
	solver_model, traj_prev = SCPS.solver_model, SCPS.traj
	robot, model = SCPP.PD.robot, SCPP.PD.model
	x_dim, u_dim, N = model.x_dim, model.u_dim, SCPP.N
	ω, Δ = SCPP.param.alg.ω_vec[end], SCPP.param.alg.Δ_vec[end]

	U = SCPV.U
	N, dt = SCPP.N, SCPP.tf_guess/SCPP.N

	cost_expr = cost_true_convexified(SCPV, traj_prev, SCPP)

	# Add penalized constraints:
	for cc_list in values(SCPC.state_trust_region_ineq)
		for cc in cc_list
			if cc.dimtype == :scalar
				cc.var_reference = @variable(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)])
				cc.con_reference = @constraint(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)],
					ω*cc.func(SCPV, traj_prev, SCPP, k, i...) - Δ <= cc.var_reference[k,i])
				@constraint(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)],
					-cc.var_reference[k,i] <= 0.)
				cost_expr += sum(cc.var_reference[k,i] for k in cc.ind_time, i in Iterators.product(cc.ind_other...))
				set_start_value.(cc.var_reference, 0.)
			elseif cc.dimtype == :array
				throw("Not implemented")
			end
		end
	end

	for cc_list in values(merge(SCPC.convex_state_ineq, SCPC.nonconvex_state_convexified_ineq))
		for cc in cc_list
			if cc.dimtype == :scalar
				cc.var_reference = @variable(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)])
				cc.con_reference = @constraint(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)],
					ω*cc.func(SCPV, traj_prev, SCPP, k, i...) <= cc.var_reference[k,i])
				@constraint(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)],
					-cc.var_reference[k,i] <= 0.)
				cost_expr += sum(cc.var_reference[k,i] for k in cc.ind_time, i in Iterators.product(cc.ind_other...))
				set_start_value.(cc.var_reference, 0.)
			elseif cc.dimtype == :array
				throw("Not implemented")
			end
		end
	end

	for cc_list in values(SCPC.convex_state_eq)
		for cc in cc_list
			if cc.dimtype == :scalar
				cc.var_reference = @variable(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...),j=1:2])
				cc.con_reference = @constraint(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...),j=1:2],
					(-1)^j*ω*cc.func(SCPV, traj_prev, SCPP, k, i...) - SCPP.param.alg.ε <= (-1)^j*cc.var_reference[k,i,j])
				@constraint(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...),j=1:2],
					-cc.var_reference[k,i,j] <= 0.)
				cost_expr += sum(cc.var_reference[k,i,j] for k in cc.ind_time, i in Iterators.product(cc.ind_other...), j in 1:2)
				set_start_value.(cc.var_reference, 0.)
			elseif cc.dimtype == :array
				throw("Not implemented")
			end
		end
	end

	@objective(solver_model, Min, cost_expr)
end

function convex_ineq_satisfied_gusto_jump(traj::Trajectory, traj_prev::Trajectory, SCPC::SCPConstraints, SCPP::SCPProblem)
  # checks for satisfaction of convex state inequalities and nonconvex->convexified state inequalities]
  for cc_list in values(merge(SCPC.convex_state_ineq, SCPC.nonconvex_state_convexified_ineq))
  	for cc in cc_list
	  	for k in cc.ind_time, i in Iterators.product(cc.ind_other...)
	  		if cc.func(traj, traj_prev, SCPP, k, i...) >= SCPP.param.alg.ε
	  			# @show cc.func, k
	  			cc.func(traj, traj_prev, SCPP, k, i...)
					return false
				end
	  	end
	  end
  end

  for cc_list in values(SCPC.convex_state_eq)
  	for cc in cc_list
	  	for k in cc.ind_time, i in Iterators.product(cc.ind_other...)
	  		constraint_value = cc.func(traj, traj_prev, SCPP, k, i...)
	  		if (constraint_value <= -SCPP.param.alg.ε) || (constraint_value >= SCPP.param.alg.ε)
	  			# @show cc.func
					return false
				end
	  	end
	  end
  end

  return true

  # TODO: Verify that this works correctly, with effect of ω and Δ
  # ω, Δ = SCPP.param.alg.ω_vec[end], SCPP.param.alg.Δ_vec[end]

  # for cc in values(merge(SCPC.convex_state_ineq, SCPC.nonconvex_state_convexified_ineq))
  # 	for k in cc.ind_time, i in Iterators.product(cc.ind_other...)
  # 		# fval = hval - rhs + slack
  # 		hval = JuMP.value(cc.con_reference[k,i])
  # 		rhs = JuMP.moi_set(JuMP.constraint_object(cc.con_reference[k,i])).upper - Δ
  # 		slack = JuMP.value(cc.var_reference[k,i])
  # 		if (hval + slack - rhs)/ω > SCPP.param.alg.ε
  # 			# @show (cc.func, k, i)
  # 			# @show (hval, rhs, slack)
  # 			# @show (SCPP.param.alg.ε, (hval + slack - rhs)/ω)
  # 			# @show (traj.U[1:3,1], norm(traj.U[1:3,1])/SCPP.PD.robot.mass)
		# 		return false
		# 	end
  # 	end
  # end
  # return true
end
