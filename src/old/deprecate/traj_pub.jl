#!/usr/bin/env julia

using RobotOS
@rosimport std_msgs.msg: Header, Float64MultiArray, MultiArrayDimension
@rosimport geometry_msgs.msg: TransformStamped, Transform, PoseStamped, Pose, Point, Point32, Vector3, Quaternion, PolygonStamped
@rosimport trajectory_msgs.msg: MultiDOFJointTrajectory, MultiDOFJointTrajectoryPoint, JointTrajectory, JointTrajectoryPoint 
rostypegen()

using std_msgs.msg
using geometry_msgs.msg
using trajectory_msgs.msg

using HDF5
using DataStructures

include("quat_functions.jl")

function main()
  init_node("body_tf_publisher")

  fn = "feasible.h5"
  Xs = h5read(fn,"Xs")

  body_transform_pub = Publisher{TransformStamped}("body_tf", queue_size=10)
  body_transform = TransformStamped()
  body_transform.header.frame_id = "map"
  body_transform.child_frame_id = "body"

  joint_trajectory = JointTrajectory()
  bay = "top_aft"
  joint_names = ["perching_arm", bay, 
                "$(@sprintf("%0.s",bay))_arm_proximal_link",
                "$(@sprintf("%0.s",bay))_arm_distal_link",
                "$(@sprintf("%0.s",bay))_gripper_left_proximal_link",
                "$(@sprintf("%0.s",bay))_gripper_left_distal_link",
                "$(@sprintf("%0.s",bay))_gripper_right_proximal_link",
                "$(@sprintf("%0.s",bay))_gripper_right_distal_link"]
  [push!(joint_trajectory.joint_names,jn) for jn in joint_names]

  loop_rate = Rate(10.0)
  idx = 1 
  while !is_shutdown()
    idx>size(Xs,2) && break

    r = vec(Xs[1:3,idx])
    q = vec(mrp2quat(Xs[4:6,idx]))

    body_transform.header.stamp = RobotOS.get_rostime()
    body_transform.transform.translation = Vector3(r[1], r[2], r[3])
    body_transform.transform.rotation = Quaternion(q[1], q[2], q[3], q[4])

    publish(body_transform_pub, body_transform)
    rossleep(loop_rate)
    idx+=1
  end
end

if !isinteractive()
  main()
end
