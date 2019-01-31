using Dierckx

function moveit_seeds(N::Int)
  fn = "trajectory_data.txt"

  joint_data = Dict()

  open(fn) do f
    line = 1
    Q = []
    key_list = []
    while !eof(f)
      x = readline(f)
      if isempty(x)
        if !isempty(key_list)
          joint_data[tuple(key_list...)] = Q
          key_list = []
          Q = []
        end
      elseif isnumber(x)
        push!(key_list, parse(Int,x))
      else
        q = [parse(Float64,ss) for ss in split(x)]
        push!(Q,q)
      end
    end
    joint_data[tuple(key_list...)] = Q;
  end

  for key in keys(joint_data)
    Qs = joint_data[key]
    n_joints = length(Qs[1])
    Qs_new = zeros(n_joints,N)

    A_x = collect(1:1:length(Qs))
    A_y = collect(linspace(1,length(Qs),N))

    for j in 1:n_joints
      q = [Q[j] for Q in Qs]
      spl = Spline1D(A_x, q, k=2)
      Qs_new[j,:] = spl(A_y)
    end
    joint_data[key] = Qs_new
  end
  return joint_data
end
