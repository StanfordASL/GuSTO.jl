using DataStructures, Graphs, MAT

function kino_rrt{T}(P::MPProblem{T},init_idx::Int=1)
  sbmp_solved = false

  # generate samples
  sample_free!(P)

  # control set
  n_controls = 12 
  Us = zeros(T,P.u_dim,n_controls)
  if P.use_2d
    warn("No controls defined")
    return
  else
    Us[1:3,1:3]   = P.robot.mass*P.robot.hard_limit_accel/sqrt(3)*eye(3)
    Us[1:3,4:6]   = -P.robot.mass*P.robot.hard_limit_accel/sqrt(3)*eye(3)
    Us[1:3,7:9]   = P.robot.mass*P.robot.hard_limit_accel/sqrt(3)*eye(3)
    Us[1:3,10:12]   = -P.robot.mass*P.robot.hard_limit_accel/sqrt(3)*eye(3)
    Us[4:6,1:3]   = P.robot.J*P.robot.hard_limit_alpha*eye(3)
    Us[4:6,4:6] = -P.robot.J*P.robot.hard_limit_alpha*eye(3)
    Us[4:6,7:9]   = -P.robot.J*P.robot.hard_limit_alpha*eye(3)
    Us[4:6,10:12] = P.robot.J*P.robot.hard_limit_alpha*eye(3)
  end

  # problem parameters
  tmin,tmax = P.t_range 
  dt_edges = tmin+(tmax-tmin)*generateHaltonSamples(1,P.n_samples)

  biases = generateHaltonSamples(1,P.n_samples)

  # create graph
  g = graph(Int[], Graphs.IEdge[], is_directed=true)
  add_vertex!(g, init_idx) # start
  nodes = Array{T}[P.init]

  tree = Int[-1]
  costs = T[0]
  times = T[0]
  ctrl_idxs = Int[0]

  tic()
  for iter = 1:P.n_samples
    x_rand = biases[iter] < P.goal_bias ? P.goal : P.V[:,iter]
    _, ind_near = findmin(map(x->dist(nodes[x],x_rand,P), g.vertices))

    # apply controls and see which guides closest to x_rand
    Xedge,min_dist,ctrl_idx = zeros(T,P.x_dim,0),T(Inf),-1
    for idx in 1:n_controls
      Xcand = simRigidEOM(nodes[ind_near],repmat(Us[:,idx],1,round(Int,dt_edges[iter]/P.dt)),P.dt,P.robot,true)
      x_new_cand = Xcand[:,end]

      dist_here = dist(x_new_cand,x_rand,P)
      if dist_here < min_dist
        Xedge = Xcand
        min_dist = dist_here 
        ctrl_idx = idx 
      end
    end

    if ctrl_idx == -1
      continue
    end

    # collision checker
    if !is_free_edge(Xedge,P)
      continue
    end

    push!(nodes, Xedge[:,end])
    push!(costs, costs[ind_near] + Us[:,ctrl_idx]'*Us[:,ctrl_idx]*dt_edges[iter])
    push!(times, dt_edges[iter])
    push!(ctrl_idxs,ctrl_idx)
    push!(tree, ind_near)

    new_id = length(nodes)
    add_vertex!(g, new_id) # start
    add_edge!(g, ind_near, new_id)

    if is_goal_pt(Xedge[:,end],P)
      print("solution found!\n")
      sbmp_solved = true
      break
    end
  end

  dim = length(nodes[1])
  Xout = zeros(T,dim, length(nodes))
  for k = 1:dim
    Xout[k,:] = map(x->nodes[x][k], 1:length(nodes))
  end

  _, c_v = findmin(map(x->dist(nodes[x],P.goal,P), g.vertices))
  path = [c_v]
  while true
    ie = in_edges(c_v,g)
    if length(ie)==0
      break
    end
    c_v = ie[1].source
    push!(path, c_v)
  end

  return reverse(path),g,tree,costs,times,Us,ctrl_idxs,Xout,sbmp_solved,toc()
end
