export PandaBot
using PandaRobot
include("../dynamics/pandadynamics.jl")

using PandaRobot, RigidBodySim, RigidBodyDynamics


mutable struct PandaBot{T<:AbstractFloat} <: Robot
  # Robot properties
  total_mass::T
  num_joints::Int
  n_motors::Int
  q_max::Vector{T}    # Joint limits
  q_min::Vector{T}
  qd_max::Vector{T}
  qd_min::Vector{T}
  qdd_max::Vector{T}
  qdd_min::Vector{T}
  qddd_max::Vector{T}
  qddd_min::Vector{T}
  tau_max::Vector{T}
  tau_min::Vector{T}
  taud_max::Vector{T}
  taud_min::Vector{T}

  # BulletPhysics struct
  btCollisionObject

  # PandaRobot struct
  pan

  q # joints # same as RigidBodyDynamics.configuration(state)
  q_dot
  q_ddot
  state # state (includes joints)

  world_frame

  EE_id::Int
  EE_link
  EE_link_frame
  EE_link_point
  EE_link_radius::T
  EE_path
  EE_point3d_inWorldframe_limit_min
  EE_point3d_inWorldframe_limit_max
end
function PandaBot{T}() where T
  total_mass = 0.0
  num_joints = 7
  n_motors = 9

  # joint space limits: https://frankaemika.github.io/docs/control_parameters.html
  q_max     = [2.8973; 1.7628; 2.8973; -0.0698; 2.8973; 3.7525; 2.8973]         # rad
  q_min     = [-2.8973; -1.7628; -2.8973; -3.0718; -2.8973; -0.0175; -2.8973]   # rad
  qd_max    = [2.1750; 2.1750; 2.1750; 2.1750; 2.6100; 2.6100; 2.6100]          # rad/s
  qdd_max   = [15;  7.5; 10; 12.5; 15; 20; 20]                                  # rad/s^2
  qddd_max  = [7500; 3750; 5000; 6250; 7500; 10000; 10000]                      # rad/s^3
  tau_max   = [87; 87; 87; 87; 12; 12; 12]                                      # N*m
  taud_max  = [1000; 1000; 1000; 1000; 1000; 1000; 1000]                        # N*m/s

  # TODO(acauligi): Instantiate Bullet representation correctly
  btCollisionObject = BulletCollision.sphere(SVector{3}(zeros(T,3)), 0.5)

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
    q0 = [0.0;0.0;0.0;0.0;0.0;pi;0.0;0.0;0.0]
    set_configuration!(state,q0)
    lims_up = RigidBodyDynamics.transform(state, EE_link_point, world_frame)
  EE_point3d_inWorldframe_limit_min, EE_point3d_inWorldframe_limit_max = -maximum(lims_up.v)*ones(3), maximum(lims_up.v)*ones(3)

  # Jacobian (joints to point)
  # EE_Jp = point_jacobian(state, EE_path, RigidBodyDynamics.transform(state, EE_link_point, world_frame))

  return PandaBot{T}(total_mass,num_joints,n_motors,
    q_max,q_min,qd_max,-qd_max,
    qdd_max,-qdd_max,qddd_max,-qddd_max,
    tau_max,-tau_max,taud_max,-taud_max,
    btCollisionObject,pan,
    q,q_dot,q_ddot,state,world_frame,
    EE_id, EE_link, EE_link_frame, EE_link_point, EE_link_radius, EE_path, 
    EE_point3d_inWorldframe_limit_min, EE_point3d_inWorldframe_limit_max)
end
PandaBot(::Type{T} = Float64; kwargs...) where {T} = PandaBot{T}(; kwargs...)


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
