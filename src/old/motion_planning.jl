#!/usr/bin/env julia

include("load_stuff.jl")
using RobotOS, MAT
using GeometryTypes

@rosimport geometry_msgs.msg: PoseStamped, Pose, Twist, Point, Quaternion, Vector3
@rosimport ff_msgs.msg: PlanActionGoal, PlanActionResult, PlanResult, ControlState, Zone, JuliaPlanRequest, JuliaPlanReply
@rosimport ff_msgs.srv: GetZones 

rostypegen()

import geometry_msgs.msg: PoseStamped, Pose, Twist, Point, Quaternion, Vector3
import ff_msgs.msg: PlanActionGoal, PlanActionResult, PlanResult, ControlState, Zone, JuliaPlanRequest, JuliaPlanReply
import ff_msgs.srv: GetZones 

function plan_cb(msg::JuliaPlanRequest,P::MPProblem)
  # See $SOURCE_PATH/communications/ff_msgs/msg/JuliaPlanRequest.msg
  # Desired state sequence
  init_pos = msg.states[1].pose.position
  init_quat = msg.states[1].pose.orientation
  P.init = [init_pos.x; init_pos.y; init_pos.z;
            zeros(3);
            quat2mrp([init_quat.x;init_quat.y;init_quat.z;init_quat.w]);
            zeros(3)]

  final_pos = msg.states[end].pose.position
  final_quat = msg.states[end].pose.orientation
  P.goal = [final_pos.x; final_pos.y; final_pos.z;
            zeros(3);
            quat2mrp([final_quat.x;final_quat.y;final_quat.z;final_quat.w]);
            zeros(3)]

  msg.faceforward        # Face-forward trajectory?
  msg.check_obstacles    # Check against obstacles?
  msg.desired_rate       # Desired rate
  msg.max_time  # Max generation time

  # Desired (max) velocity
  P.robot.hard_limit_vel = msg.desired_vel

  # Desired (max) accel
  P.robot.hard_limit_accel = msg.desired_accel

  # Desired (max) omega
  P.robot.hard_limit_omega = msg.desired_omega

  # Desired (max) alpha
  P.robot.hard_limit_alpha = msg.desired_alpha

  println("Plan received in Julia node!")
  P.received_plan=true
end

function parse_zones!{T}(P::MPProblem{T})
  # See $SOURCE_PATH/communications/ff_msgs/msg/Zone.ms
  const srvcall = ServiceProxy("/mob/get_zones", GetZones)
  wait_for_service("/mob/get_zones")
  zones_srv = srvcall(ff_msgs.srv.GetZonesRequest())
  
  P.world.keepin_zones = Vector{HyperRectangle}(0)
  P.world.keepout_zones = Vector{HyperRectangle}(0)
  for (k,zone) in enumerate(zones_srv.zones)
    lo = [zone.min.x, zone.min.y, zone.min.z]
    hi = [zone.max.x, zone.max.y, zone.max.z]
    new_zone = rectangle(lo,hi)

    # types- KEEPOUT:0, KEEPIN:1, CLUTTER:2
    push!(
      zone.type==UInt8(1) ? P.world.keepin_zones : P.world.keepout_zones, 
      HyperRectangle(Vec3f0(lo),Vec3f0(hi-lo))
    )
  end

  matwrite("ISSzones.mat", Dict(
    "keepin_zones" => P.world.keepin_zones,
    "keepout_zones" => P.world.keepout_zones
  ))

  update_aabb!(P.world)
  nzones = length(P.world.keepout_zones)+length(P.world.keepin_zones)
  println("Total zones received: $nzones")
end

function package_controls{T}(Xs::Matrix{T},Us::Matrix{T},Ts::Vector{T})
  reply = JuliaPlanReply()

  for k in 1:size(Xs,2)
    state = ControlState()
    state.when = Time(Ts[k])
    x = Xs[:,k]

    state.pose.position.x = x[1]
    state.pose.position.y = x[2]
    state.pose.position.z = x[3]
    state.twist.linear = Vector3(x[4],x[5],x[6])

    q = mrp2quat(x[7:9])
    state.pose.orientation = Quaternion(q[1],q[2],q[3],q[4])
    
    state.twist.angular = Vector3(x[10],x[11],x[12])

    k!=size(Xs,2) && continue
    if k != 1
      f,m  = quat_rotate(q,Us[1:3,k-1]), Us[4:6,k-1]
      state.accel.linear = Vector3(f[1],f[2],f[3]) 
      state.accel.angular = Vector3(m[1],m[2],m[3])
    end
    push!(reply.segment, state)
  end
  return reply
end

function dynProp{T}(P::MPProblem{T},X0::Vector{T},Us::Matrix{T},Ts::Vector{T},dt::T)
  # X0: initial 12-dim vector
  # Us: control applied as ZOH
  # Ts: time stamps for controls applied
  Xs = repmat(X0,1,length(Ts))
  x0 = copy(X0)

  for k in 1:length(Ts)-1
    U = Us[:,k]
    Us = repmat(Us,1,Int(round((Ts[k+1]-Ts[k])/dt)))
    Xout = simRigidEOM(x0,Us,dt,P.robot)
    Xs[:,k+1] = Xout[:,end]
    x0 = Xout[:,end]
  end
  return Xs 
end

function dummy_control{T}(P::MPProblem{T})
  reply = JuliaPlanReply()
  N = 10
  Ts = collect(linspace(0,9,N))

  U = P.use_2d ? 
    0.25*P.robot.hard_limit_accel/sqrt(2)*[0;0;0;1;1;0.] :
    0.25*P.robot.hard_limit_accel/sqrt(3)*[0;0;0;1;1;1.] 
  Us = repmat(U,1,N-1)

  Xs = dynProp(P,P.init,Us,Ts,P.dt)

  for k in 1:size(Xs,2)
    state = ControlState()
    state.when = Ts[k] 
    x = Xs[:,k]

    state.pose.position.x = x[1]
    state.pose.position.y = x[2]
    state.pose.position.z = x[3]
    state.twist.linear = Vector3(x[4],x[5],x[6])

    q = mrp2quat(x[7:9])
    state.pose.orientation = Quaternion(q[1],q[2],q[3],q[4])
    
    state.twist.angular = Vector3(x[10],x[11],x[12])

    k!=size(Xs,2) && continue
    if k != 1
      f,m  = quat_rotate(q,Us[1:3,k-1]), Us[4:6,k-1]
      state.accel.linear = Vector3(f[1],f[2],f[3]) 
      state.accel.angular = Vector3(m[1],m[2],m[3])
    end
    push!(reply.segment, state)
  end

  for k in 1:size(Xs,2)
    println("$(Xs[7,k]) $(Xs[8,k]) $(Xs[9,k]) $(Xs[10,k])")
  end

  return reply 
end

function main()
  init_node("julia_planner_fmt",anonymous=true)
  P = MPProblem()

  sub = Subscriber{JuliaPlanRequest}("julia_plan_request", plan_cb, (P,), queue_size=10)
  pub = Publisher{JuliaPlanReply}("julia_plan_reply", queue_size=10)
  parse_zones!(P)

  loop_rate = Rate(5.0)
  while !is_shutdown() && !P.received_plan
    println("Waiting...")
    rossleep(loop_rate)
  end

  # params defined in $SOURCE_PATH/astrobee/config/mobility/planner_fmt.config
  P.use_2d = get_param("/cfg/planner_fmt/two_d")
  P.robot.r = get_param("/cfg/planner_fmt/robot_radius")

  # path,g,tree,costs,times,controls,ctrl_idxs,Xout,sbmp_solved,krrt_time = kino_rrt(P)
  # Dt = P.segment_dt
  # Tf = sum(times[path])
  # Xs = P.init 
  # Us = zeros(P.u_dim,0)
  # for (k,path_idx) in enumerate(path[2:end])
  #   X0 = Xs[:,end]
  #   N = round(Int,times[path_idx]/P.dt)
  #   newUs = repmat(controls[:,ctrl_idxs[path_idx]],1,N)
  #   newXs = simRigidEOM(X0,newUs,P.dt,P.robot,true)[:,2:end]
  #   Xs = [Xs newXs]
  #   Us = [Us newUs]
  # end
  # Ts = collect(linspace(0,sum(times[path]),size(Xs,2)))

  # println("Tf")
  # for k in 1:size(Xs,2)
  #   # println("$(Xs[7,k]) $(Xs[8,k]) $(Xs[9,k])")
  #   println("$() $() $() $() $() $()")
  # end

  # # downsample to 2Hz
  # f = 2

  planning_reply = package_controls(Xs,Us,Ts)

  publish(pub,package_controls(Xs,Us,Ts))
  println("publishing")
end

if !isinteractive()
  main()
end
