export PandaBot
using PandaRobot
include("../dynamics/pandadynamics.jl")

using PandaRobot, RigidBodySim, RigidBodyDynamics

mutable struct panda_bubble
  parent_id::Int
  local_pose::Vector
  radius::Number
end

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

  # collision bubbles
  bubble_array

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

  bubble_array = panda_bubbles()
  objs = Vector{BulletCollision.BulletCollisionObjectPtr}(0)
  for bubble in bubble_array
    push!(objs,
      BulletCollision.sphere(SVector{3}(zeros(T,3)), bubble.radius))
  end
  btCollisionObject = BulletCollision.compound_collision_object(objs)

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

  return PandaBot{T}(total_mass,num_joints,n_motors,
    q_max,q_min,qd_max,-qd_max,
    qdd_max,-qdd_max,qddd_max,-qddd_max,
    tau_max,-tau_max,taud_max,-taud_max,
    bubble_array,btCollisionObject,pan,
    q,q_dot,q_ddot,state,world_frame,
    EE_id, EE_link, EE_link_frame, EE_link_point, EE_link_radius, EE_path)
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
# Collision info 
##################
function panda_bubbles()
  # Information for spheres along kinematic chain stored here
  arr = Vector{panda_bubble}(0)

  # Link 1
  push!(arr, panda_bubble(3,[0.;0.;-0.12],0.08))
  push!(arr, panda_bubble(3,[-0.0;-0.05;-0.02],0.08))

  # Link 2
  push!(arr, panda_bubble(4,[0.;-0.02;0.06],0.07))
  push!(arr, panda_bubble(4,[0.;-0.08;0.03],0.07))
  push!(arr, panda_bubble(4,[0.;-0.14;0.01],0.07))

  # Link 3
  push!(arr, panda_bubble(5,[0.;0.;-0.08],0.05))
  push!(arr, panda_bubble(5,[0.04;0.04;-0.03],0.05))
  push!(arr, panda_bubble(5,[0.1;0.07;0.02],0.05))

  # Link 4
  push!(arr, panda_bubble(6,[0.;0.;0.06],0.06))
  push!(arr, panda_bubble(6,[-0.06;0.02;0.04],0.06))
  push!(arr, panda_bubble(6,[-0.07;0.08;0.01],0.06))

  # Link 5
  push!(arr, panda_bubble(7,[0.;0.01;-0.2],0.06))
  push!(arr, panda_bubble(7,[0.;0.05;-0.13],0.06))
  push!(arr, panda_bubble(7,[0.;0.09;-0.06],0.06))
 
  # Link 6
  push!(arr, panda_bubble(8,[0.0;0.0;-0.08],0.05))
  push!(arr, panda_bubble(8,[-0.02;0.0;0.02],0.03))
  push!(arr, panda_bubble(8,[0.01;0.0;0.02],0.03))
  
  # Link 7
  push!(arr, panda_bubble(9,[0.0;0.0;0.1],0.05))
  push!(arr, panda_bubble(9,[0.0;0.0;0.15],0.04))
  push!(arr, panda_bubble(9,[0.04;0.05;0.15],0.04))
  push!(arr, panda_bubble(9,[-0.04;-0.05;0.15],0.04))
  return arr
end
