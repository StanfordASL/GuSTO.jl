using Primes
using Convex, SCS
using MAT

export generateHaltonSamples, yershova_quaternions

function generateHaltonSamples(dim::Int, n::Int, sf::Array=[])
  primeNums = primes(1000)

  if length(sf) != dim
    # warn("Scale factors for Halton samples being set to 1 by default!")
    sf = ones(dim)
  end

  samples = zeros(dim, n)
  
  for (idx, base_val) in enumerate(primeNums[1:dim]) 
    samples[idx, :] = sf[idx] > 0 ? sf[idx]*generateHaltonSequence(n, base_val) : -2*sf[idx]*(generateHaltonSequence(n, base_val) - 0.5)
  end
  return samples
end

function generateHaltonSequence(n::Int, base::Int) 
  halton_sequence = zeros(n, 1)

  for idx in range(1, 1, n)
    halton_sequence[idx] = localHaltonSingleNumber(idx, Float64(base))
  end
  return halton_sequence
end

function localHaltonSingleNumber(n::Int, b::Float64)
  n0 = n
  hn = 0
  f = 1/b

  while(n0 > 0)
    n1 = floor(n0/b)
    r = n0 - b*n1
    hn = hn + f*r
    f = f/b
    n0 = n1
  end
  return hn
end

function yershova_quaternions()
  # http://rotations.mitchell-lab.org/
  # Yershova, Anna, et al. "Generating uniform incremental grids on SO (3) using the Hopf fibration." 
  # The International journal of robotics research 29.7 (2010): 801-812.
  return readdlm(joinpath(ENV["GUSTO"], "quat_data.txt"))
end

function uniform_random_quaternion(n::Int, canonicalize::Bool=false)
  # Algorithm 2 from "Effective Sampling and Distance Metrics for 3D Rigid Body Path Planning" by James Kuffner (2004)
  q_matrix = zeros(4,n)

  for idx in 1:n
    s = rand()
    sigma1 = sqrt(1-s)
    sigma2 = sqrt(s)

    theta1 = 2*pi*rand()
    theta2 = 2*pi*rand()

    x = sin(theta1) * sigma1
    y = cos(theta1) * sigma1
    z = sin(theta2) * sigma2
    w = cos(theta2) * sigma2
   
    temp_q = [x y z w]
    temp_q = canonicalize ? sign(w)*temp_q./norm(temp_q) : temp_q./norm(temp_q)

    q_matrix[:,idx] = temp_q
  end
  return q_matrix
end

function uniform_random_euler(n::Int)
  # Algorithm 1 from "Effective Sampling and Distance Metrics for 3D Rigid Body Path Planning" by James Kuffner (2004)
  euler_angles = zeros(3,n)
  idx = 1
  
  while idx <= n
    roll = 2 * pi * rand() - pi
    pitch = acos(1 - 2* rand()) + pi/2

    if rand() < 0.5
      if pitch < pi
        pitch = pitch+pi
      elseif pitch == pi
        continue
      else
        pitch = pitch-pi;
      end
    end

    yaw = 2*pi*rand() - pi
    
    euler_angles[:, idx] = [roll; pitch; yaw]
    idx+=1
  end

  return euler_angles
end

# function control_samples(n_samples::Int, rb::astrobee)
#   # writes matrix to file with 6 x n_samples control vectors
#   # each control vector multiplied by robot mass and inertia i.e. [F;M]
#   # f defined in body frame, [F;M] used in inertial frame
#   # but norm of G*f is preserved independent of orientation
# 
#   n_samples_max = 5*n_samples
#   Us = [generateHaltonSamples(3, n_samples_max, -rb.hard_limit_accel/sqrt(3)*[1;1;1]); 
#         generateHaltonSamples(3, n_samples_max, -rb.hard_limit_alpha/sqrt(3)*[1;1;1])]
# 
#   controls = zeros(Float64,6,0)
#   fs = zeros(Float64,rb.n_thrusters,0)
# 
#   f = Variable(rb.n_thrusters,1)
#   prob = minimize(norm(rb.G*f,1))
#   prob.constraints += 0 <= f
#   prob.constraints += norm2(rb.G[1:3,:]*f) <= rb.hard_limit_accel
#   prob.constraints += norm2(rb.G[4:6,:]*f) <= rb.hard_limit_alpha
# 
#   Convex.solve!(prob, SCSSolver(verbose=false))
# 
#   f_idx = 1
#   @time for U_idx = 1:size(Us,2)
#     prob.constraints +=  Us[:,U_idx] == rb.G*f
#     Convex.solve!(prob, SCSSolver(verbose=false))
#     
#     if prob.status == :Optimal
#       U = [rb.mass*ones(3); diag(rb.J)].*rb.G*f.value
#       controls = [controls U]
#       fs = [fs f.value]
#       f_idx+=1
# 
#       if f_idx > n_samples
#         break
#       end
#     end
#     prob.constraints = prob.constraints[1:end-1]
#   end
# 
#   file = matopen("control_samples.mat","w")
#   write(file, "Us", controls)
#   close(file)
# end
