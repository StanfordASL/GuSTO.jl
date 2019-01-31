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
