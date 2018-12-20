export PandaBot
using PandaRobot, RigidBodySim, RigidBodyDynamics

mutable struct PandaBot{T<:AbstractFloat} <: Robot
  q_max
  q_min
  qd_max
  qdd_max
  qddd_max
  tau_max
  tau_min
  taud_max
  taud_min
  mechanism
  btCollisionObject
end
function PandaBot{T}() where T
  # joint space limits: https://frankaemika.github.io/docs/control_parameters.html
  q_max     = [2.8973; 1.7628; 2.8973; -0.0698; 2.8973; 3.7525; 2.8973]         # rad
  q_min     = [-2.8973; -1.7628; -2.8973; -3.0718; -2.8973; -0.0175; -2.8973]   # rad
  qd_max    = [2.1750; 2.1750; 2.1750; 2.1750; 2.6100; 2.6100; 2.6100]          # rad/s
  qdd_max   = [15;  7.5; 10; 12.5; 15; 20; 20]                                  # rad/s^2
  qddd_max  = [7500; 3750; 5000; 6250; 7500; 10000; 10000]                      # rad/s^3
  tau_max   = [87; 87; 87; 87; 12; 12; 12]                                      # N*m
  tau_min   = -tau_max                                                          # N*m
  taud_max  = [1000; 1000; 1000; 1000; 1000; 1000; 1000]                        # N*m/s
  taud_min  = -taud_max                                                         # N*m/s

  mechanism = Panda()
  btCollisionObject = BulletCollision.sphere(SVector{3}(zeros(T,3)), 0.5)

  return PandaBot{T}(q_max,q_min,qd_max,qdd_max,qddd_max,tau_max,tau_min,taud_max,taud_min,
    mechanism,btCollisionObject)
end
PandaBot(::Type{T} = Float64; kwargs...) where {T} = PandaBot{T}(; kwargs...)

# function misc()
#   pd = Panda()
#   state = RigidBodyDynamics.MechanismState(pd.mechanism)
#   world = RigidBodyDynamics.root_frame(pd.mechanism)
#   
#   pd_root_body = root_body(pd.mechanism)
#   pd_hand_body = findbody(pd.mechanism, "panda_hand")
#   point = Point3D(default_frame(pd_hand_body), 0.,0.,0.)
#   p = RigidBodyDynamics.path(pd.mechanism, pd_root_body, pd_hand_body)
# 
#   # forward kinematics
#   after_panda_hand_joint = RigidBodyDynamics.Spatial.transform(state,point,world)
#   after_panda_hand_joint.v
# 
#   # Jacobian
#   Jp = point_jacobian(state, p, after_panda_hand_joint)
#   Jp.J
# end
