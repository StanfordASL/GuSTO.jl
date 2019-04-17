mutable struct Goal{T<:GoalType}
	params::T

	t_guess 				# Goal time guess
	ind_coordinates # Coordinate indices on which goal is defined
	k_timestep 			# Assigned timestep number
	t_final  				# Final calculated goal time
end

function Goal(params, t_guess, ind_coordinates)
	Goal(params, t_guess, ind_coordinates, missing, missing)
end

function Goal(params, t_guess, model::M) where M <: DynamicsModel
	Goal(params, t_guess, 1:model.x_dim, missing, missing)
end

function assign_timesteps!(goal_set::GoalSet, N, tf_guess)
	for (t_guess, goal) in goal_set.goals
		goal.k_timestep = Int(fld(N*tf_guess, N*t_guess))
	end
end

function get_first_goal_at_time(goal_set::GoalSet, t)
	deref_value((goal_set.goals, searchsortedfirst(goal_set.goals, t)))
end

function add_goal!(goal_set::GoalSet, goal::Goal{T}) where T
	insert!(goal_set.goals, goal.t_guess, goal)
end

mutable struct PointGoal <: GoalType
	point
end

mutable struct BoxGoal <: GoalType
	lower_bound
	upper_bound
end

mutable struct BallGoal <: GoalType
	center
	radius
end

center(goal::PointGoal) = goal.point
center(goal::BoxGoal) = 1/2*(goal.upper_bound + goal.lower_bound)
center(goal::BallGoal) = goal.center