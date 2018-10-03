include("load_common.jl")
using MAT

using AstrobeeRobot

using RigidBodySim
using RigidBodyDynamics

using MeshCat
using MeshCatMechanisms
import GeometryTypes: HyperRectangle, Vec, HomogenousMesh
import ColorTypes: RGBA

#=  ----------------------------
	---- DISTANCE FUNCTIONS ----
 	---------------------------- =#  
include("bullet_collision.jl")


# # compute distances to all objects and save results in a table
# if(!isdefined(:DIST2OBS_MAT))
#     NX = 131 # number of points along each dimension
#     X_MIN, X_MAX = x0-(xf-x0)/4, xf+(xf-x0)/4
#     DX = (X_MAX[1]-X_MIN[1])
#     dx = DX/(NX-1)
#     println("The distances table will have a precision of $(dx/2) meters")
#     grid_precision = dx/2
  # SEP11_2018: OUTDATED FUNCTION, MIGHT NOT WORK AS INTENDED
  # SEP11_2018: OUTDATED FUNCTION, MIGHT NOT WORK AS INTENDED
  # SEP11_2018: OUTDATED FUNCTION, MIGHT NOT WORK AS INTENDED
  #  function create_table_distances(NX_::Any, X_MIN_::Any, DX_::Any)
  #      last_println = 0
  #      println("[create_table_distancesxid] Starting to compute the distances table...")
  #      for xid = 1:NX_
  #          if last_println>=NX_/10
  #              println("[create_table_distancesxid] Done at $(100*xid/NX_)%")
  #              last_println = 0
  #          else
  #              last_println = last_println + 1
  #          end
  #
  #          x = X_MIN_[1] + (xid-1)*DX_/(NX_-1)
  #          for yid = 1:NX_
  #              y = X_MIN_[2] + (yid-1)*DX_/(NX_-1)
  #              for zid = 1:NX_
  #                  z = X_MIN_[3] + (zid-1)*DX_/(NX_-1)
  #                  pt_coords_3d = [x,y,z]
  #                  min_distance = Inf
  #                  for obs in P.world.btenvironment_keepout.convex_env_components
  #                      pt = BT.convex_hull([SVector{3}(pt_coords_3d)]) # represent a point using convex hull
  #                      distance_to_obs = BT.distance(pt,obs)[1]# zero if inside obstacle
  #                      if min_distance > distance_to_obs
  #                          min_distance = distance_to_obs
  #                      end
  #                  end
  #                  DIST2OBS_MAT[xid,yid,zid] = min_distance
  #              end
  #          end
  #      end   
  #      println("[create_table_distancesxid] Computed distance table")
  #      return DIST2OBS_MAT
  #  end
#     # table with all signed euclidean distances to the closest obstacle for all points
#     DIST2OBS_MAT = zeros(NX,NX,NX)
#     DIST2OBS_MAT = create_table_distances();
# end

# for a specific set of obstacles
function create_table_distances(NX_::Any, X_MIN_::Any, DX_::Any, bt_obstacles_::Any)
    last_println = 0
    println("[create_table_distancesxid] Starting to compute the distances table...")
    for xid = 1:NX_
        if last_println>=NX_/10
            println("[create_table_distancesxid] Done at $(100*xid/NX_)%")
            last_println = 0
        else
            last_println = last_println + 1
        end

        x = X_MIN_[1] + (xid-1)*DX_/(NX_-1)
        for yid = 1:NX_
            y = X_MIN_[2] + (yid-1)*DX_/(NX_-1)
            for zid = 1:NX_
                z = X_MIN_[3] + (zid-1)*DX_/(NX_-1)
                pt_coords_3d = [x,y,z]
                min_distance = Inf
                for obs in bt_obstacles_.convex_env_components
                    pt = BT.convex_hull([SVector{3}(pt_coords_3d)]) # represent a point using convex hull
                    distance_to_obs = BT.distance(pt,obs)[1]# zero if inside obstacle
                    if min_distance > distance_to_obs
                        min_distance = distance_to_obs
                    end
                end
                DIST2OBS_MAT[xid,yid,zid] = min_distance
            end
        end
    end   
    println("[create_table_distancesxid] Computed distance table")
    return DIST2OBS_MAT
end


function get_dist2obs{T}(pt_3d_::Vector{T}, DIST2OBS_MAT_::Array{T,3}, NX_::Any, X_MIN_::Any, DX_::Any)
    idx_in_tab = zeros(3);
    for dim in 1:3
        idx_in_tab[dim] = (pt_3d_[dim]-X_MIN_[dim]) / DX_
        if idx_in_tab[dim] < 0
#             println("not enough")
            idx_in_tab[dim] = 1
        elseif idx_in_tab[dim] > 1
#             println("too much")
            idx_in_tab[dim] = NX_
        else
            idx_in_tab[dim] = 1+idx_in_tab[dim]*(NX_-1)
        end
    end
    idx_in_tab = round.(Int,idx_in_tab)
    # x = X_MIN_[1] + (idx_in_tab[1]-1)*DX_/(NX_-1)
    # y = X_MIN_[2] + (idx_in_tab[2]-1)*DX_/(NX_-1)
    # z = X_MIN_[3] + (idx_in_tab[3]-1)*DX_/(NX_-1)
    #println("index in tab: $idx_in_tab located at [$x $y $z]")
    return DIST2OBS_MAT_[idx_in_tab[1],idx_in_tab[2],idx_in_tab[3]]
end

#=  ------------------------------------
	---- VELOCITY/TORQUES FUNCTIONS ----
 	------------------------------------ =#  
function get_finite_differencing_acceleration_matrix(N_::Any, dt_::Any)
    Aacc_ = zeros(N_,N_) # to take into account the constant part!
    Aacc_ = Aacc_ - 30.0/12.0*diagm(ones(N_))
    Aacc_ = Aacc_ + 16.0/12.0*diagm(ones(N_-1),1) + 16.0/12.0*diagm(ones(N_-1),-1)
    Aacc_ = Aacc_ - 1/12.0*diagm(ones(N_-2),2) - 1/12.0*diagm(ones(N_-2),-2)
    Aacc_ = Aacc_/(dt_^2)
    return Aacc_
end
 	# Computes the velocity magnitudes for each 3d point in 'points'
function compute_velocities_mag{T}(points::Matrix{T})
	p = size(points,1)
	N = size(points,2)
#     diff_mat_vel = zeros(N,N)
#     diff_mat_vel = diagm(ones(N)) - diagm(ones(N-1),-1); diff_mat_vel[1,1] = 0;
#     diff_mat_vel = diff_mat_vel / dt
    diff_mat_vel = zeros(N+1,N+1)
    diff_mat_vel = diagm(ones(N+1)) - diagm(ones(N+1-1),-1);
    diff_mat_vel = diff_mat_vel / dt
    
    velocities = similar(points);
    pts_extended = zeros(p,N+1)
    pts_extended = hcat(points[:,1],points[:,:])
    for dim=1:p
#         velocities[dim,:] = (diff_mat_vel * points[dim,:])'
        velocities[dim,:] = (diff_mat_vel * pts_extended[dim,:])'[2:end]
    end

    vel_mag = zeros(N)
    for i=1:N
        vel_mag[i] = sqrt(sum(velocities[:,i].^2))
    end
    
    return vel_mag
end
# Compute all velocities
# velocities_mag = zeros(N)
# velocities_mag = compute_velocities_mag(Q0)

# Computes the torques required for each motor (cost)
function compute_torques{T}(Q_::Matrix{T}, NB_MOTORS_::Any, body_mass_::Any, body_inertia_::Matrix{T}) # configurations
    # TODO CHANGE THAT, REPLACE WITH INVERSE DYNAMICS ALGORITHMS, see eq (15)
	p = size(Q_,1)
	N = size(Q_,2)
	Aacc = get_finite_differencing_acceleration_matrix(N+2*2, dt)
        
    # put inverse dynamics here
    # put inverse dynamics here
    # put inverse dynamics here
    Q_accels = zeros(p,N)
    for dim=1:p
        # Assumes initial points are the same as 1st element of Q, and similarly for the end
#         Q_accels[dim,:] = A * Q[dim, :]
        Q_accels[dim,:] = (Aacc * [Q_[dim,1];Q_[dim,1];Q_[dim,:];Q_[dim,end];Q_[dim,end]])[3:end-2]
    end
    torques = zeros(NB_MOTORS_,N) # torques of each motor for each timestep
    for nidx=1:N
        for motor_i = 1:NB_MOTORS_
            if motor_i <= 3
#                 torques[motor_i,:] = body_mass_*sqrt(sum(Q_accels[:,nidx]).^2)
                torques[motor_i,nidx] = body_mass_*sqrt((Q_accels[motor_i,nidx]).^2)
            elseif motor_i > 3 && motor_i <= p && p > 3
#                 torques[motor_i,nidx] = body_inertia_[*sqrt(sum(Q_accels[motor_i,nidx]).^2)
#                torques[motor_i,nidx] = body_inertia_[motor_i-3,motor_i-3]*sqrt((Q_accels[motor_i,nidx]).^2)
                torques[motor_i,nidx] = body_inertia_[1,1]*sqrt((Q_accels[motor_i,nidx]).^2)
            end
        end
    end
    # put inverse dynamics here
    # put inverse dynamics here
    # put inverse dynamics here
        
    return torques
end
# Compute all velocities
# torques = zeros(NB_MOTORS,N)
# torques = compute_torques(Q0)


#=  ----------------------------
	----   COST FUNCTIONS   ----
 	---------------------------- =#  
function get_total_jointsAccels_cost{T}(Q_::Matrix{T})
    Jacc = 0
#     Jacc1 = 0
    for dim = 1:size(Q_,1)
        # Assumes initial points are the same as 1st element of Q, and similarly for the end
#         Jacc1 = Jacc1 + Q_[dim,:]' * R * Q_[dim,:]
        Qextended = [Q_[dim,1];Q_[dim,1];Q_[dim,:];Q_[dim,end];Q_[dim,end]]
        Jacc = Jacc + Q_[dim,:]' * ((Racc * Qextended)[3:end-2])
    end
    
    return Jacc
end
# get_total_jointsAccels_cost(Q0)
function get_total_obstacles_cost{T}(Q_::Matrix{T}, eps_::Any, body_radius_::Any, DIST2OBS_MAT_::Array{T,3}, NX_::Any, X_MIN_::Any, DX_::Any)
#     N = T
    
    velocities_mag_ = zeros(size(Q,2))
    velocities_mag_ = compute_velocities_mag(Q_)
    
    Jobs_ = 0
    for nidx = 1:size(Q_,2)
        # put forward kinematics here
        # put forward kinematics here
        qconf_ = Q_[1:3, nidx] 
        # put forward kinematics here
        # put forward kinematics here

        dist_ = get_dist2obs(qconf_, DIST2OBS_MAT_, NX_, X_MIN_, DX_)
        Jobs_ = Jobs_ + max(eps_+body_radius_-dist_, 0)*velocities_mag_[nidx]
    end
    
    return Jobs_
end
function is_trajectory_collision_less{T}(Q_::Matrix{T}, eps_::Any, body_radius_::Any, DIST2OBS_MAT_::Array{T,3}, NX_::Any, X_MIN_::Any, DX_::Any)
    is_not_collision = true
    for nidx = 1:size(Q_,2)
        # put forward kinematics here
        # put forward kinematics here
        qconf_ = Q_[1:3, nidx] 
        # put forward kinematics here
        # put forward kinematics here

        dist_ = get_dist2obs(qconf_, DIST2OBS_MAT_, NX_, X_MIN_, DX_)
        if (eps_+body_radius_-dist_ > 0)
            is_not_collision = false
            return is_not_collision
        end
    end
    
    return is_not_collision
end
function get_total_torques_cost{T}(Q_::Matrix{T}, NB_MOTORS_::Any, weight_torques_::Any, body_mass_::Any, body_inertia_::Matrix{T})
	N = size(Q_,2)
    
    velocities_mag_ = zeros(N)
    velocities_mag_ = compute_velocities_mag(Q_)
    torques_ = zeros(NB_MOTORS_,N)
    torques_ = compute_torques(Q_, NB_MOTORS_, body_mass_, body_inertia_)
    
    Jtorques_ = 0
    for nidx = 1:N
        Jtorques_ = Jtorques_ + weight_torques_*sum(abs.(torques_[:,nidx]))
    end
    
    return Jtorques_
end
function get_total_refDeviation_cost{T}(Q_::Matrix{T}, Q0_::Matrix{T}, weight_ref_::Any)
	N = size(Q_,2)
    
    Jref_ = 0
    for nidx = 1:N
        qconf_ = Q_[:, nidx] 
        Jref_ = Jref_ + weight_ref_*sqrt(sum((qconf_-Q0_[:,nidx]).^2))
    end
    
    return Jref_
end

# Cost function. Returns the cost for a given configuration qconf
# and the velocity of the corresponding qconf identified by id_vel (precomputed in velocities_mag)
function cost{T}(qconf_::Vector{T}, Q0_::Matrix{T}, id_vel_::Any, nidx_::Any, eps_obs_::Any, body_radius_::Any, 
				 velocities_mag_::Vector{T}, torques_::Matrix{T}, weight_torques_::Any, weight_ref_::Any, 
				 DIST2OBS_MAT_::Array{T,3}, NX_::Any, X_MIN_::Any, DX_::Any)

	state.q = qconf_
	EE_point_inWorldFrame = RigidBodyDynamics.transform(state, EE_link_point, world)


    dist_ = get_dist2obs(EE_point_inWorldFrame.v, DIST2OBS_MAT_, NX_, X_MIN_, DX_)
    
    # obstacle cost (eq 13)
#     println(dist)
    cost_obs_ = max(eps_obs_+body_radius_-dist_, 0)*velocities_mag_[id_vel_]
    
    # joint torques costs
    cost_torques_ = weight_torques_*sum(abs.(torques_[:,nidx_]))
#     cost_torques_ = 0
    
    # distance to reference cost
    cost_qconf_ref_ = weight_ref_*sqrt(sum((qconf_-Q0_[:,nidx_]).^2))
#     cost_qconf_ref_ = 0
    
    Jm_ = 1e-5 + cost_obs_ + cost_torques_ + cost_qconf_ref_
    return Jm_
end

# Cost function. Returns the cost for a given configuration qconf
# and the velocity of the corresponding qconf identified by id_vel (precomputed in velocities_mag)
function cost_println{T}(qconf_::Vector{T}, Q0_::Matrix{T}, id_vel_::Any, nidx_::Any, eps_obs_::Any, body_radius_::Any, 
						 velocities_mag_::Vector{T}, torques_::Matrix{T}, weight_torques_::Any, weight_ref_::Any, 
						 DIST2OBS_MAT_::Array{T,3}, NX_::Any, X_MIN_::Any, DX_::Any)
    dist_ = get_dist2obs(qconf_, DIST2OBS_MAT_, NX_, X_MIN_, DX_)
    
    # obstacle cost (eq 13)
#     println(dist)
    cost_obs_ = max(eps_obs_+body_radius_-dist_, 0)*velocities_mag_[id_vel_]
    
    # joint torques costs
    cost_torques_ = weight_torques_*sum(abs.(torques_[:,nidx_]))
#     cost_torques_ = 0
    
    # distance to reference cost
    cost_qconf_ref_ = weight_ref_*sqrt(sum((qconf_-Q0_[:,nidx_]).^2))
#     cost_qconf_ref_ = 0
    
    Jm_ = 1e-5 + cost_obs_ + cost_torques_ + cost_qconf_ref_
    println("Contribution of costs: ")
    println("cost_obs:       $(cost_obs/Jm)%")
    println("cost_torques:   $(cost_torques/Jm)%")
    println("cost_qconf_ref: $(cost_qconf_ref/Jm)%")
    return Jm_
end