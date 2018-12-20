export Plane
  
mutable struct Plane{T<:AbstractFloat} <: Robot
  mass::T 
  g::T
  ρ::T
  Area::T
  Cd0::T
  Kd::T
  α_0::T
  v_hat::T    # z_dot \in [-0.1,0.1]*v_hat
  v_min::T    # v \in [3,10] m/s
  v_max::T
  γ_max::T    # gamma \in [-π/6,π/6]
  ϕ_max::T    # phi \in [-π/4,π/4]
  acc_lim::T  # u_a \in [-8,8] m/s
  ϕ_d_lim::T  # ϕ_d \in [-120,120] deg/s
  α_d_lim::T  # α_d \in [-60,60] deg/s
  btCollisionObject
end

function Plane{T}() where T
  # X: (x,y,z,psi,v,gamma,phi, alpha)
  # pos- xyz, heading- psi, speed- v, flight path angle- gamma, roll- phi, angle of attack- alpha,
  # U: (u_acc, u_ϕ_d, u_α_d)

  mass = 1.
  g = 9.81
  ρ = 1.225
  Area = 0.7
  Cd0 = 0.015
  Kd = 0.025
  α_0 = 5*π/180
  γ_max = π/6
  ϕ_max = π/4
  v_min = 3.    # min speed
  v_max = 10.   # max speed

  C_con = π*ρ*Area/mass
  v_hat = sqrt(g/(C_con*α_0))
  acc_nom = (ρ*Area)*(v_hat^2)*(Cd0+4*(π^2)*Kd*(α_0^2))

  acc_lim = 10*acc_nom      # linear acceleration lim
  ϕ_d_lim = (π/3)*2      # ϕ_dot lim
  α_d_lim = 2*γ_max # α_dot lim

  btCollisionObject = BulletCollision.sphere(SVector{3}(zeros(3)), 1.)

  return Plane{T}(mass,g,ρ,Area,Cd0,Kd,α_0,v_hat,v_min,v_max,γ_max,ϕ_max,acc_lim,ϕ_d_lim,α_d_lim,btCollisionObject)
end
Plane(::Type{T} = Float64; kwargs...) where {T} = Plane{T}(; kwargs...)

function configuration(x, rb::Plane)
  # convert state to (x,y,z,yaw,pitch,roll)
  pitch = rb.α_0 - x[8] - x[6]
  return [x[1:4]; pitch; x[7]]
end

function body2inertial(x,rb)
  # See 2.2 in Small Unmanned Aircraft: Theory and Practice by Beard & McLain

  X = configuration(x,rb) # X: (x,y,z,yaw,pitch,roll)

  R_z = [cos(yaw) -sin(yaw) 0;
         sin(yaw) cos(yaw) 0;
         0 0 1]
     
  pitch = -pitch  # since not using NED
  
  R_y = [cos(pitch) 0 sin(pitch);
         0 1 0;
        -sin(pitch) 0 cos(pitch)]
    
  R_x = [1 0 0;
         0 cos(roll) -sin(roll);
         0 sin(roll) cos(roll)]
    
  R = R_z*R_y*R_x   # R: body->inertial
end

function inertial2body(x,rb)
  return body2inertial(x,rb)' # R: inertial->body
end

function plane_geometry()
  α_0 = 5*π/180
  λ = 0.1
  
  # Vertices
  body_V =
      [10. -0.5 -1.
       10.  0.5 -1.
       10. -0.5  0.  
       10.  0.5  0.
        5. -0.5  1.  
        5.  0.5  1.
      -10. -0.5  1.  
      -10.  0.5  1.
      -10. -0.5  0. 
      -10.  0.5  0.
      -7.5 -0.5 -1. 
      -7.5  0.5 -1.]*λ
  body = FlexibleConvexHull([Vector(body_V[k,:]) for k in 1:size(body_V,1)])
  
  wing_V = [ 3.   0.   3*sin(α_0)+.1
             3.   0.   3*sin(α_0)-.1
             0. -10.  .1
             0. -10. -.1
             0.  10.  .1
             0.  10. -.1
            -2. -10.  -2*sin(α_0)+.1
            -2. -10.  -2*sin(α_0)-.1
            -2.  10.  -2*sin(α_0)+.1
            -2.  10.  -2*sin(α_0)-.1]*λ
  wing = FlexibleConvexHull([Vector(wing_V[k,:]) for k in 1:size(wing_V,1)])
      
  stab_V = [ -8.  0. .5
             -8.  0. .1
             -9. -4. .5  
             -9. -4. .1
             -9.  4. .5  
             -9.  4. .1
             -10. -4. .5
             -10. -4. .1
             -10.  4. .5
             -10.  4. .1]*λ
  stab = FlexibleConvexHull([Vector(stab_V[k,:]) for k in 1:size(stab_V,1)])
  
  tail_V =  [-7. -.2 1.
             -6. .2 1.
              -10. -.2 1.
              -10. .2 1.
              -10. -.2 4.
              -10. .2 4.
              -9. -.2 4.
              -9. .2 4.]*λ
  tail = FlexibleConvexHull([Vector(tail_V[k,:]) for k in 1:size(tail_V,1)])
  
  box_V = [-1 -1 -0.1
          1 -1 -0.1
          1 1 -0.1
          -1 1 -0.1
          -1 -1 0.4
          1 -1 0.4
          1  1 0.4
          -1 1 0.4]
  box = FlexibleConvexHull([Vector(box_V[k,:]) for k in 1:size(box_V,1)])
  
  return body, wing, stab, tail, box
end
