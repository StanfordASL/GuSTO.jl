export rrt
using DataStructures, Graphs, MAT

include("tpbvp.jl")

function dist{T}(v1::Vector{T},v2::Vector{T},P::MPProblem)
  dot = v1[7:10]'*v2[7:10]
  if (abs(dot) > 1.)
    dot = sign(dot)
  end
  return norm(v1[1:3]-v2[1:3]) + P.slew_weighting*acos(dot)
end

function rrt!{T}(P::MPProblem{T})
  tic()
  sbmp_solved = false
  init_idx = 1
  init = P.init[1:10]

  goal_bias = 0.05
  eps = 0.5 

  biases = generateHaltonSamples(1,P.n_samples)

  # create graph
  g = graph(Int[], Graphs.IEdge[], is_directed=true)
  add_vertex!(g, init_idx) # start
  c_v = init_idx # closest vertex
  nodes = Array{T}[P.init]

  tree = Int[-1]
  costs = T[0]

  tic()
  for k = 1:P.n_samples
    x_rand = biases[k] < goal_bias ? P.goal : P.V[:,k]

    _, ind_near = findmin(map(x->dist(nodes[x],x_rand,P), g.vertices))

    eps_here = min(1., eps/dist(nodes[ind_near],x_rand,P))
    x_near = [nodes[ind_near][1:6]+ eps_here*(x_rand[1:6]-nodes[ind_near][1:6]); 
              quat_interp(nodes[ind_near][7:10], x_rand[7:10], eps_here)]

    # collision checker
    # if !BT.is_free_motion()
    if false
      continue
    end

    push!(nodes, x_near)
    new_id = length(nodes)

    push!(costs, costs[ind_near] + dist(nodes[ind_near],x_rand,P))
    push!(tree, ind_near)

    add_vertex!(g, new_id) # start
    add_edge!(g, ind_near, new_id)

    if is_goal_pt(x_near,P)
      print("solution found!\n")
      sbmp_solved = true
      break
    end
  end

  Xout = zeros(length(init), length(nodes))
  for k = 1:length(init)
    Xout[k,:] = map(x->nodes[x][k], 1:length(nodes))
  end

  _, c_v = findmin(map(x->norm(nodes[x][1:3]-P.goal[1:3]), g.vertices))
  path = [c_v]
  while true
    ie = in_edges(c_v,g)
    if length(ie)==0
      break
    end
    c_v = ie[1].source
    push!(path, c_v)
  end

  P.solution = MPSolution(reverse(path),tree,T[],costs,Xout,
                            sbmp_solved,toc(),false,T(0),[])
  return
end

function plotter(path,Xout)
  for k = 2:length(path)
    idx1,idx2 = path[k-1], path[k]
    PyPlot.plot3D([Xout[1,idx1], Xout[1,idx2]],
          [Xout[2,idx1], Xout[2,idx2]],
          [Xout[3,idx1], Xout[3,idx2]])
  end
end
