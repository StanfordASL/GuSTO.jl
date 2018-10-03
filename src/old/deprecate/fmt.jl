export fmtstar!
using DataStructures

function fmtstar!{T}(P::MPProblem{T}; init_idx=1, checkpts=false)
  tic()

  times         = P.tpbvp_data.times
  adj_costs     = P.tpbvp_data.adj_costs
  nn_go_edges   = P.tpbvp_data.nn_go_edges
  nn_come_edges = P.tpbvp_data.nn_come_edges

  if false 
    warn("Initial state is infeasible!")
    P.solution = MPSolution(Int[],Int[],T[],T[],[],
                            false,toc(),false,T(0),[])
    return
  end

  if checkpts # TODO(acauligi): collision checking
    F = trues(P.n_samples)
    # for i in 1:P.n_samples
    #   F[i] = BT.is_free_state(P.V[:,i])
    # end
  end

  A = zeros(Int,P.n_samples)    # parent
  W = trues(P.n_samples)        # unvisited
  H = falses(P.n_samples)       # frontier
  C = zeros(T,P.n_samples)      # costs
  W[init_idx] = false
  H[init_idx] = true
  HHeap = binary_minheap(Pair{T,Int})
  push!(HHeap, Pair(T(0), init_idx))

  cost,z = pop!(HHeap)   # i.e. z = init_idx

  ctr = 1
  while !is_goal_pt(z,P)
    H_new = Int[]

    for x in nn_go_edges[P.k_nearest*(z-1)+1:P.k_nearest*z]
      !W[x] && continue # node has been visited 
      (x<1||x>P.n_samples) && continue
      checkpts && !F[x] && continue

      c_min,y_min = C[z]+adj_costs[P.n_samples*(z-1)+x],z
      for y in nn_come_edges[P.k_nearest*(x-1)+1:P.k_nearest*x]
        !H[y] && continue # node not on frontier
        (y<1||y>P.n_samples) && continue

        if (C[y]+adj_costs[P.n_samples*(y-1)+x] < c_min)
          c_min = C[y]+adj_costs[P.n_samples*(y-1)+x]
          y_min = y
        end
      end

      if true # TODO(acauligi): collision checking
        A[x] = y_min
        C[x] = c_min
        push!(HHeap, Pair(c_min, x))
        push!(H_new, x)
        W[x] = false
      end
    end
    H[H_new] = true
    H[z] = false
    if !isempty(HHeap)
      cost,z = pop!(HHeap)
    else
      break
    end
  end

  sbmp_solved = is_goal_pt(z,P)

  sol = [z]
  costs = [C[z]]
  while sol[1] != 1
    unshift!(sol, A[sol[1]])
    if sol[1] == -1
      unshift!(costs, T(0))
      break
    end
    unshift!(costs, C[sol[1]])
  end

  soln_times = Vector{T}()
  soln_costs = Vector{T}()
  for (k,src) in enumerate(sol[1:end-1])
    dst = sol[k+1]
    idx = P.n_samples*(src-1)+dst 
    push!(soln_times,P.tpbvp_data.times[idx])  
    push!(soln_costs,P.tpbvp_data.adj_costs[idx])  
  end

  P.solution = MPSolution(sol,[],soln_times,soln_costs,P.V[:,sol],
                sbmp_solved,toc(),false,T(0),[])
end
