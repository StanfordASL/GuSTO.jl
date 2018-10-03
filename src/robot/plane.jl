export Plane
  
mutable struct Plane{T<:AbstractFloat} <: Robot
  mass::T 
  g::T
  rho::T
  Area::T
  Cd0::T
  Kd::T
  alpha_0::T
  v_hat::T        # z_dot \in [-0.1,0.1]*v_hat
  v_min::T        # v \in [3,10] m/s
  v_max::T
  gamma_max::T    # gamma \in [-pi/6,pi/6]
  phi_max::T      # phi \in [-pi/4,pi/4]
  acc_lim::T      # u_a \in [-8,8] m/s
  phi_d_lim::T    # phi_d \in [-120,120] deg/s
  alpha_d_lim::T  # alpha_d \in [-60,60] deg/s
  btCollisionObject
end

function Plane{T}() where T
  # X: (x,y,z,psi,v,gamma,phi, alpha)
  # pos- xyz, heading- psi, speed- v, flight path angle- gamma, roll- phi, angle of attack- alpha,
  # U: (u_acc, u_phi_d, u_alpha_d)

  mass = 1.
  g = 9.81
  rho = 1.225
  Area = 0.7
  Cd0 = 0.015
  Kd = 0.025
  alpha_0 = 5*pi/180
  gamma_max = pi/6
  phi_max = pi/4
  v_min = 3.    # min speed
  v_max = 10.   # max speed

  C_con = pi*rho*Area/mass
  v_hat = sqrt(g/(C_con*alpha_0))
  acc_nom = (rho*Area)*(v_hat^2)*(Cd0+4*(pi^2)*Kd*(alpha_0^2))

  acc_lim = 10*acc_nom      # linear acceleration lim
  phi_d_lim = (pi/3)*2      # phi_dot lim
  alpha_d_lim = 2*gamma_max # alpha_dot lim

  btCollisionObject = BulletCollision.sphere(SVector{3}(zeros(3)), 1.)

  return Plane{T}(mass,g,rho,Area,Cd0,Kd,alpha_0,v_hat,v_min,v_max,gamma_max,phi_max,acc_lim,phi_d_lim,alpha_d_lim,btCollisionObject)
end
Plane(::Type{T} = Float64; kwargs...) where {T} = Plane{T}(; kwargs...)

function configuration(x, rb::Plane)
  # convert state to (x,y,z,yaw,pitch,roll)
  pitch = rb.alpha_0 - x[8] - x[6]
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
  alpha_0 = 5*pi/180
  lambda = 0.1
  
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
      -7.5  0.5 -1.]*lambda
  body = FlexibleConvexHull([Vector(body_V[k,:]) for k in 1:size(body_V,1)])
  
  wing_V = [ 3.   0.   3*sin(alpha_0)+.1
             3.   0.   3*sin(alpha_0)-.1
             0. -10.  .1
             0. -10. -.1
             0.  10.  .1
             0.  10. -.1
            -2. -10.  -2*sin(alpha_0)+.1
            -2. -10.  -2*sin(alpha_0)-.1
            -2.  10.  -2*sin(alpha_0)+.1
            -2.  10.  -2*sin(alpha_0)-.1]*lambda
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
             -10.  4. .1]*lambda
  stab = FlexibleConvexHull([Vector(stab_V[k,:]) for k in 1:size(stab_V,1)])
  
  tail_V =  [-7. -.2 1.
             -6. .2 1.
              -10. -.2 1.
              -10. .2 1.
              -10. -.2 4.
              -10. .2 4.
              -9. -.2 4.
              -9. .2 4.]*lambda
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
