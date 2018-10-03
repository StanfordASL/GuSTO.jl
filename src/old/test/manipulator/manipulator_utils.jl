# Manipulator utils

using MeshCat
using MeshCatMechanisms

include("../../TrajectoryOptimization/robot/pandabot.jl")


# Inverse Kinematics
# from http://nbviewer.jupyter.org/github/JuliaRobotics/RigidBodyDynamics.jl/blob/master/notebooks/Jacobian%20IK%20and%20Control.ipynb
function jacobian_transpose_ik!(state::MechanismState, # current state (used as initial guess)
                               body::RigidBody,  # body we are using
                               point::Point3D,   # point we want to set to desired
                               desired::Point3D, # desired position of the point
                               α::T=0.1,
                               iterations::Int=100) where T
    mechanism = state.mechanism
    world = root_frame(mechanism)
    
    # Compute the joint path from world to our target body
    p = RigidBodyDynamics.path(mechanism, root_body(mechanism), body)
    # Allocate the point jacobian (we'll update this in-place later)
    Jp = point_jacobian(state, p, RigidBodyDynamics.transform(state, point, world))
    
    q = copy(RigidBodyDynamics.configuration(state))
    
    for i in 1:iterations
        # Update the position of the point
        point_in_world = RigidBodyDynamics.transform(state, point, world)
        # Update the point's jacobian
        point_jacobian!(Jp, state, p, point_in_world)
        # Compute an update in joint coordinates using the jacobian transpose
        Δq = α * Array(Jp)' * (RigidBodyDynamics.transform(state, desired, world) - point_in_world).v
        # Apply the update
        q .= RigidBodyDynamics.configuration(state) .+ Δq
        set_configuration!(state, q)
        #println(i)
    end
    state
end
# Inverse Kinematics
# from http://nbviewer.jupyter.org/github/JuliaRobotics/RigidBodyDynamics.jl/blob/master/notebooks/Jacobian%20IK%20and%20Control.ipynb
function inverse_kinematics!(manipulator::PandaBot,#state::MechanismState, # current state (used as initial guess)
                             body::RigidBody,  # body we are using
                             point::Point3D,   # point we want to set to desired
                             desired::Point3D, # desired position of the point
                             α::T=0.1,
                             iterations::Int=200) where T
    desired_in_world = RigidBodyDynamics.transform(manipulator.state, desired, manipulator.world_frame)
    if maximum(manipulator.EE_point3d_inWorldframe_limit_min .> desired_in_world.v) ||
       maximum(manipulator.EE_point3d_inWorldframe_limit_max .< desired_in_world.v)
       println("[inverse_kinematics!] WARNING: Point not reachable.")
    end
    ## TODO: RETURN ERROR IF TOO BIG ERROR (IF UNREACHABLE)
    ## TODO: RETURN ERROR IF TOO BIG ERROR (IF UNREACHABLE)
    ## TODO: RETURN ERROR IF TOO BIG ERROR (IF UNREACHABLE)
    mechanism = manipulator.state.mechanism
    world = root_frame(mechanism)
    
    # Compute the joint path from world to our target body
    p = RigidBodyDynamics.path(mechanism, root_body(mechanism), body)
    # Allocate the point jacobian (we'll update this in-place later)
    Jp = point_jacobian(manipulator.state, p, RigidBodyDynamics.transform(manipulator.state, point, world))
    
    q = copy(RigidBodyDynamics.configuration(manipulator.state))
    
    for i in 1:iterations
        # Update the position of the point
        point_in_world = RigidBodyDynamics.transform(manipulator.state, point, world)
        # Update the point's jacobian
        point_jacobian!(Jp, manipulator.state, p, point_in_world)
        # Compute an update in joint coordinates using the jacobian transpose
        Δq = α * Array(Jp)' * (RigidBodyDynamics.transform(manipulator.state, desired, world) - point_in_world).v
        # Apply the update
        q .= RigidBodyDynamics.configuration(manipulator.state) .+ Δq
        set_configuration!(manipulator.state, q)
        #println(i)
    end

    update_q_robot(manipulator)
    ## TODO: RETURN ERROR IF TOO BIG ERROR (IF UNREACHABLE)
    ## TODO: RETURN ERROR IF TOO BIG ERROR (IF UNREACHABLE)
    ## TODO: RETURN ERROR IF TOO BIG ERROR (IF UNREACHABLE)

    manipulator.state
end

