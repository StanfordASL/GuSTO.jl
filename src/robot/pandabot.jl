export PandaBot
using PandaRobot
include("../dynamics/pandadynamics.jl")

using RigidBodySim
using RigidBodyDynamics


mutable struct PandaBot{T<:AbstractFloat} <: Robot
  # Robot Properties
  total_mass::T
  n_motors::Int
  # J::Matrix{T}
  # Jinv::Matrix{T}
  # n_thrusters::Int
  r::T
  # hard_limit_vel::T
  # hard_limit_accel::T 
  # hard_limit_omega::T
  # hard_limit_alpha::T
  btCollisionObject

  # PandaRobot structure
  pan#::Panda <: PandaRobot

  # Robot state
  q # joints # same as RigidBodyDynamics.configuration(state)
  q_dot
  q_ddot
  state # state (includes joints)

  # Frame
  world_frame

  # End effector
  EE_id::Int
  EE_link
  EE_link_frame
  EE_link_point
  EE_link_radius::T
  EE_path
  EE_point3d_inWorldframe_limit_min
  EE_point3d_inWorldframe_limit_max
  # EE_Jp
end
function PandaBot{T}() where T
  # Robot Properties
  total_mass = 0.0
  n_motors = 9
  # s = 0.5*0.305   # each side of cube is 30.5cm
  r = sqrt(3)*0.5*0.305   # inflate to sphere
  # hard_limit_vel = 0.5 
  # hard_limit_accel = 0.1 
  # hard_limit_omega = 45*π/180 
  # hard_limit_alpha = 50*π/180 
  # J = 0.1083*Eye(3) 
  # Jinv = inv(J)
  btCollisionObject = BT.sphere(SVector{3}(zeros(T,3)), r) # Update this

  # PandaRobot structure
  # Note that it includes the definition of all the links and their relative definition
  pan = Panda()

  # Robot state
  q = zeros(n_motors)
  q_dot = zeros(n_motors)
  q_ddot = zeros(n_motors)
  state = RigidBodyDynamics.MechanismState(pan.mechanism)
    set_configuration!(state,q)
    set_velocity!(state,q_dot)


  # Frame
  world_frame = RigidBodyDynamics.root_frame(pan.mechanism)

  # End effector (hand)
  EE_id = 11
  EE_link = RigidBodyDynamics.bodies(pan.mechanism)[EE_id]
  EE_link_frame = RigidBodyDynamics.default_frame(EE_link)
  EE_link_point = RigidBodyDynamics.Point3D(EE_link_frame,0.0,0.0,0.1)
  EE_link_radius = 0.02
  EE_path = RigidBodyDynamics.path(pan.mechanism, root_body(pan.mechanism),EE_link)
     # Get limits
    q0 = [0.0;0.0;0.0;0.0;0.0;π;0.0;0.0;0.0]
    set_configuration!(state,q0)
    lims_up = RigidBodyDynamics.transform(state, EE_link_point, world_frame)
  EE_point3d_inWorldframe_limit_min, EE_point3d_inWorldframe_limit_max = -maximum(lims_up.v)*ones(3), maximum(lims_up.v)*ones(3)
  # Jacobian (joints to point)
  # EE_Jp = point_jacobian(state, EE_path, RigidBodyDynamics.transform(state, EE_link_point, world_frame))  )

  # new PandaBot instance
  return PandaBot{T}(total_mass,n_motors,r,btCollisionObject,pan,q,q_dot,q_ddot,state,world_frame,
    EE_id, EE_link, EE_link_frame, EE_link_point, EE_link_radius, EE_path, 
    EE_point3d_inWorldframe_limit_min, EE_point3d_inWorldframe_limit_max)
  # return PandaBot{T}(total_mass,n_motors,pan,q,state)
end
PandaBot(::Type{T} = Float64; kwargs...) where {T} = PandaBot{T}(; kwargs...)


##################
# Initializations
##################
# Careful, x_init, x_goal are in configuration space!
function init_traj_straightline(TOP::TrajectoryOptimizationProblem{PandaBot{T}, PandaDyn, E}) where {T,E}
  model, x_init, x_goal = TOP.PD.model, TOP.PD.x_init, TOP.PD.x_goal
  x_dim, u_dim, N, tf_guess = model.x_dim, model.u_dim, TOP.N, TOP.tf_guess
  N = TOP.N

  X = hcat(linspace(x_init, x_goal, N)...)
  U = zeros(u_dim, N)
  Trajectory(X, U, tf_guess)
end


##################
# Access Values
##################
function get_EE_point_inWorldFrame(manipulator::PandaBot)
  EE_point_inWorldFrame = RigidBodyDynamics.transform(manipulator.state, manipulator.EE_link_point, manipulator.world_frame)
  return EE_point_inWorldFrame
end
function get_mass_matrix(manipulator::PandaBot)
  mass_matrix(manipulator.state)
end
function get_torques(manipulator::PandaBot)
  inverse_dynamics(manipulator.state, manipulator.q_dot)
end

##################
# Update Values
##################
function update_state_robot(manipulator::PandaBot)
      set_configuration!(manipulator.state, manipulator.q)
      # set_configuration!(mvis, configuration(state))
end
function update_q_robot(manipulator::PandaBot)
      manipulator.q = manipulator.state.q
      # set_configuration!(mvis, configuration(state))
end
