export BlankEnv

mutable struct BlankEnv <: Environment 
  worldAABBmin::Vector
  worldAABBmax::Vector
  keepin_zones::Vector
  keepout_zones::Vector
  obstacle_set::Vector
end

function BlankEnv()
  worldAABBmin = -1000*ones(3)
  worldAABBmax = 1000*ones(3)
  keepin_zones = Vector{GeometryTypes.GeometryPrimitive}(0)
  keepout_zones = Vector{GeometryTypes.GeometryPrimitive}(0)
  obstacle_set = Vector{GeometryTypes.GeometryPrimitive}(0)
  BlankEnv(worldAABBmin, worldAABBmax, keepin_zones, keepout_zones, obstacle_set)
end



function add_sphere_random!(env::BlankEnv, 
                            p_init, p_goal, 
                            radius::T=0.05) where T<:AbstractFloat
  # dont take points 0.3% close to either p_init or p_goal, only in the middle
  center = zeros(3)
  for i=1:3
      center[i] = p_init[i] + (p_goal[i]-p_init[i])*(0.4*rand()+0.3)
  end
  center = Vec3f0(center[1], center[2], center[3])

  push!(env.obstacle_set, HyperSphere(Point3f0(center), Float32(radius)))
  # ----------
end

function sample_blank_env_with_n_random_spheres(Nb_envs::Int, Nb_obs::Int, 
                                                p_init, p_goal, 
                                                radius::T=0.05) where T<:AbstractFloat
  envs = []

  for i = 1:Nb_envs
    env   = BlankEnv();

    for j = 1:Nb_obs
      add_sphere_random!(env, p_init, p_goal)
    end

    push!(envs, env)
  end

  return envs
end