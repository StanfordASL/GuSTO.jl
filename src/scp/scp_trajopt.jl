export solve_trajopt_jump!

mutable struct SCPParam_TrajOpt <: SCPParamSpecial
  mu0         # initial penalty coefficient
  s0          # initial trust region size
  c           # step acceptance parameter
  τ_plus      # trust region expansion factor
  τ_minus     # trust region shrinkage factor
  k           # penalty scaling factor
  ftol        # convergence threshold for merit
  xtol        # convergence threhsold for x
  ctol        # constraint satisfaction threshold
  max_penalty_iteration
  max_convex_iteration
  max_trust_iteration

  ρ_vec::Vector   # trust region ratios
  mu_vec::Vector  # penalty coefficient
  s_vec::Vector   # trust region size
  xtol_vec::Vector
  ftol_vec::Vector
  ctol_vec::Vector
end

function SCPParam_TrajOpt(mu0, s0, c, τ_plus, τ_minus, k, ftol, xtol, ctol,max_penalty_iteration,max_convex_iteration,max_trust_iteration)
	SCPParam_TrajOpt(mu0, s0, c, τ_plus, τ_minus, k, ftol, xtol, ctol, 
    max_penalty_iteration,max_convex_iteration,max_trust_iteration,[0.],[mu0], [s0], [0.], [0.], [0.])
end

#########
# JuMP 
#########
function solve_trajopt_jump!(SCPS::SCPSolution, SCPP::SCPProblem, solver="Ipopt", max_iter=125, force=false; kwarg...)
	# Solves a sequential convex programming problem
	# Inputs:
	#		SCPS - SCP solution structure, including an initial trajectory
	#		SCPP - SCP problem definition structure
	#   max_iter - Set to enforce a maximum number of iterations (0 = unlimited)
	#   force - Set true to force further refinement iterations despite convergence

	N = SCPP.N
	param, model = SCPP.param, SCPP.PD.model

	param.alg = SCPParam_TrajOpt(SCPP.PD.model)
  mu0,s0 = param.alg.mu0, param.alg.s0
  ρ_vec,mu_vec,s_vec = param.alg.ρ_vec, param.alg.mu_vec, param.alg.s_vec
  τ_plus,τ_minus = param.alg.τ_plus, param.alg.τ_minus
  c = param.alg.c
  ftol,xtol,ctol = param.alg.ftol, param.alg.xtol, param.alg.ctol
  ftol_vec, xtol_vec, ctol_vec = param.alg.ftol_vec, param.alg.xtol_vec, param.alg.ctol_vec

  max_penalty_iteration,max_convex_iteration,max_trust_iteration = 
    param.alg.max_penalty_iteration, param.alg.max_convex_iteration, param.alg.max_trust_iteration

  iter_cap = SCPS.iterations + max_iter
	SCPV = SCPVariables{JuMP.VariableRef, Array{JuMP.VariableRef}}()
	SCPS.SCPC = SCPConstraints(SCPP)
	SCPC = SCPS.SCPC

  initialize_model_params!(SCPP, SCPS.traj)
	push!(SCPS.J_true, cost_true(SCPS.traj, SCPS.traj, SCPP))
	param.obstacle_toggle_distance = model.clearance + 1. # TODO(ambyld): Generalize this, workspace-dependent
	
  constraints_satisfied = false
  xtol_satisfied = false

  old_penalty_traj, old_convex_traj = SCPS.traj, SCPS.traj

  for penalty_iteration in 1:max_penalty_iteration
    constraints_satisfied && break
    copy!(old_penalty_traj,SCPS.traj)

    for convex_iteration in 1:max_convex_iteration
      copy!(old_convex_traj,SCPS.traj)

      constraints_satisfied && break
      if xtol_satisfied
        xtol_satisfied = false
        break
      end
      for trust_iteration in 1:max_trust_iteration
        time_start = time_ns()

        if solver == "Mosek"
          SCPS.solver_model = Model(with_optimizer(Mosek.Optimizer; kwarg...))
        elseif solver == "Gurobi"
          SCPS.solver_model = Model(with_optimizer(Gurobi.Optimizer; kwarg...))
        elseif solver == "Ipopt"
          SCPS.solver_model = Model(with_optimizer(Ipopt.Optimizer; kwarg...))
        else 	# Default solver
          SCPS.solver_model = Model(with_optimizer(Ipopt.Optimizer; kwarg...))
        end

        # Set up, solve problem
        update_model_params!(SCPP, SCPS.traj)
        add_variables_jump!(SCPS, SCPV, SCPP)
        add_constraints_trajopt_jump!(SCPS, SCPV, SCPC, SCPP)   # TODO(acauligi)
        add_objective_trajopt_jump!(SCPS, SCPV, SCPC, SCPP)     # TODO(acauligi)

        set_start_value.(SCPV.X, SCPS.traj.X)
        set_start_value.(SCPV.U, SCPS.traj.U)
        set_start_value(SCPV.Tf, SCPP.tf_guess)

        JuMP.optimize!(SCPS.solver_model)

		    push!(SCPS.solver_status, JuMP.termination_status(SCPS.solver_model))
        if SCPS.solver_status[end] ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
          warn("TrajOpt iteration failed to find an optimal solution")
          push!(SCPS.iter_elapsed_times, (time_ns() - time_start)/10^9) 
        end

        # Recover solution
		    new_traj = Trajectory(JuMP.value.(SCPV.X), JuMP.value.(SCPV.U), JuMP.value.(SCPV.Tf))
        push!(xtol_vec, evaluate_xtol(new_traj, old_convex_traj, SCPP))
        push!(SCPS.convergence_measure, xtol_vec[end])
		    push!(SCPS.J_full, JuMP.objective_value(SCPS.solver_model))

		    SCPS.dual = get_dual_jump(SCPC, SCPP)

        # grow or shrink trust region 
        push!(ρ_vec, trust_region_ratio_trajopt(new_traj, old_convex_traj, SCPP))
        if ρ_vec[end] > c
          push!(s_vec, τ_plus*s_vec[end])
        else
          push!(s_vec, τ_minus*s_vec[end])
        end
        
        copy!(SCPS.traj, new_traj)
        iter_elapsed_time = (time_ns() - time_start)/10^9
        push!(SCPS.J_true, cost_true(new_traj, SCPS.traj, SCPP))
        push!(SCPS.iter_elapsed_times, iter_elapsed_time)
        SCPS.total_time += iter_elapsed_time
        SCPS.iterations += 1
        if s_vec[end] < xtol
          xtol_satisfied = true
          break
        end
      end
      
      push!(ftol_vec, evaluate_ftol(SCPS.traj, old_convex_traj, SCPP))
      push!(xtol_vec, evaluate_xtol(SCPS.traj, old_convex_traj, SCPP))
      if ftol_vec[end] < ftol || xtol[end] < xtol
        constraints_satisfied = true
        break
      end
    end

    if evaluate_ctol(SCPS.traj, old_penalty_traj, SCPP, SCPC) < ctol
  	 constraints_satisfied = true
     SCPS.converged = true
     break
  	else 
      # increase penalty coefficient 
      push!(mu_vec, mu_vec[end]*SCPP.param.alg.k)
    end
  end
end

function add_constraints_trajopt_jump!(SCPS::SCPSolution, SCPV::SCPVariables, SCPC::SCPConstraints, SCPP::SCPProblem)
  s = SCPP.param.alg.s_vec[end]
	solver_model, traj_prev = SCPS.solver_model, SCPS.traj

	update_model_params!(SCPP, traj_prev)

  for cc_list in values(SCPC.state_trust_region_ineq)
    for cc in cc_list
  		if cc.dimtype == :scalar
  			cc.con_reference = @constraint(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)], cc.func(SCPV, traj_prev, SCPP, k, i...) - s <= 0)
  		elseif cc.dimtype == :array
  			cc.con_reference = @constraint(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)], cc.func(SCPV, traj_prev, SCPP, k, i...) - s .<= 0)
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

  for cc_list in values(SCPC.convex_state_eq)
    for cc in cc_list
  		if cc.dimtype == :scalar
  			cc.con_reference = @constraint(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)], cc.func(SCPV, traj_prev, SCPP, k, i...) == 0)
  		elseif cc.dimtype == :array
  			cc.con_reference = @constraint(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)], cc.func(SCPV, traj_prev, SCPP, k, i...) .== 0)
  		end
    end
	end

  # TODO: Store constraint reference to this
	if !SCPP.param.fixed_final_time
		@constraint(solver_model, SCPV.Tf >= 0.1)
	end
end


function add_objective_trajopt_jump!(SCPS::SCPSolution, SCPV::SCPVariables, SCPC::SCPConstraints, SCPP::SCPProblem)
	solver_model, traj_prev = SCPS.solver_model, SCPS.traj
	robot, model = SCPP.PD.robot, SCPP.PD.model
	x_dim, u_dim, N = model.x_dim, model.u_dim, SCPP.N

  mu = SCPP.param.alg.mu_vec[end]

	cost_expr = cost_true_convexified(SCPV, traj_prev, SCPP)

  # Inequality constraints
  for cc_list in values(merge(SCPC.convex_state_ineq,SCPC.nonconvex_state_convexified_ineq,SCPC.convex_control_ineq))
    for cc in cc_list
  		if cc.dimtype == :scalar
  			cc.var_reference = @variable(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)])
  			cc.con_reference = @constraint(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)],
  				mu*cc.func(SCPV, traj_prev, SCPP, k, i...) <= cc.var_reference[k,i])
  			@constraint(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)],
  				-cc.var_reference[k,i] <= 0.)
  			cost_expr += sum(cc.var_reference[k,i] for k in cc.ind_time, i in Iterators.product(cc.ind_other...))
  			set_start_value.(cc.var_reference, 0.)
  		elseif cc.dimtype == :array
  			throw(":array type inequality constraint penalty not implemented")      # TODO(acauligi)
  		end
    end
	end

  # Equality constraints
  for cc_list in values(SCPC.nonconvex_state_convexified_eq)
    for cc in cc_list
  		if cc.dimtype == :scalar
  			cc.var_reference = @variable(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)])
  			cc.con_reference = @constraint(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)],
  				mu*cc.func(SCPV, traj_prev, SCPP, k, i...) <= cc.var_reference[k,i])
        @constraint(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)],
  				-mu*cc.func(SCPV, traj_prev, SCPP, k, i...) <= cc.var_reference[k,i])         # TODO(acauligi): Need to track this second constraint satisfaction somehow
  			@constraint(solver_model, [k=cc.ind_time,i=Iterators.product(cc.ind_other...)],
  				-cc.var_reference[k,i] <= 0.)
  			cost_expr += sum(cc.var_reference[k,i] for k in cc.ind_time, i in Iterators.product(cc.ind_other...))
  			set_start_value.(cc.var_reference, 0.)
  		elseif cc.dimtype == :array
  			throw(":array type equality constraint penalty not implemented")      # TODO(acauligi)
  		end
    end
	end

  # Dynamics constraints
  for cc_list in values(SCPC.dynamics)
    for cc in cc_list
      if cc.dimtype == :scalar
        throw(":scalar type dynamics constraint penalty not implemented")      # TODO(acauligi)
      elseif cc.dimtype == :array
        # TODO(ambyld): Redesign this
        for k in cc.ind_time
          var_reference1 = @variable(solver_model, [j=1:x_dim])
          var_reference2 = @variable(solver_model, [j=1:x_dim])
          @constraint(solver_model, mu*cc.func(SCPV, traj_prev, SCPP, k) - var_reference1 .<= 0)
          @constraint(solver_model, var_reference2 - mu*cc.func(SCPV, traj_prev, SCPP, k) .<= 0)
          cost_expr += sum(var_reference1[j] + var_reference2[j] for j in 1:x_dim)
          set_start_value.(var_reference1, 0.)
          set_start_value.(var_reference2, 0.)
        end
      end
    end
  end


	@objective(solver_model, Min, cost_expr) 
end

function evaluate_xtol(traj::Trajectory, traj_prev::Trajectory, SCPP::SCPProblem)
  return convergence_metric(traj,traj_prev,SCPP)
end

function evaluate_ftol(traj::Trajectory, traj_prev::Trajectory, SCPP::SCPProblem)
  abs(cost_true(traj, traj, SCPP)-cost_true(traj_prev, traj_prev, SCPP))/abs(cost_true(traj, traj, SCPP))
end

function evaluate_ctol(traj::Trajectory, traj_prev::Trajectory, SCPP::SCPProblem, SCPC::SCPConstraints)
	JNum, JDen = 0, 0
  JtestNum, JtestDen = 0, 0

  fprev = SCPC.convex_state_ineq[1].func
  for cc_list in values(merge(SCPC.convex_state_ineq, SCPC.nonconvex_state_ineq, SCPC.convex_state_boundary_condition_ineq, SCPC.nonconvex_state_boundary_condition_ineq, SCPC.convex_state_eq, SCPC.dynamics, SCPC.nonconvex_state_eq))
    for cc in cc_list
      fnext = cc.func
      for k in cc.ind_time, i in Iterators.product(cc.ind_other...)
        if fnext == fprev	# Evaluate and update maximum value while considering same class of constraints
          JtestNum = max(JtestNum, norm(fnext(traj, traj, SCPP, k, i...) -fnext(traj_prev, traj_prev, SCPP, k, i...)))
          JtestDen = max(JtestDen, norm(fnext(traj, traj, SCPP, k, i...)))
        else		# Entering new class of constraints
          JNum += JtestNum	# Add max from previous class of constraints to cost
          JDen += JtestDen
          fprev = fnext	# Update function handle representing current class of constraints
          JtestNum = norm(fnext(traj, traj, SCPP, k, i...) - fnext(traj_prev, traj_prev, SCPP, k, i...)) 	# Evaluate value of first member of new class of constraints
          JtestDen = norm(fnext(traj, traj, SCPP, k, i...))
        end
      end
    end
  end
	return JNum/JDen
end
