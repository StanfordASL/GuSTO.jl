using NLsolve, DifferentialEquations, Convex, SCS, PyPlot

function dubins_dyn(X,U,dt) 
  x,y,th = X
  w = U

  f = [V*cos(th);
      V*sin(th);
      w]
  Ak = [0 0 -V*sin(th);
        0 0 V*cos(th);
        0 0 0]
  Bk = [0;0;1]
  return f,Ak,Bk
end

function dubins(dx,x,p,t)
  x,y,th,px,py,pth = x
  u = k*pth/2.0
  dx[1] = V*cos(th)
  dx[2] = V*sin(th)
  dx[3] = k*u
  dx[4] = 0.
  dx[5] = 0.
  dx[6] = px*V*sin(th) - py*V*cos(th)
end

function shooting!(Xout,lambda0) 
  X0 = [x0;lambda0]
  tspan = (0.0,Tf)
  prob = ODEProblem(dubins,X0,tspan)
  sol = solve(prob)
	
  Xout = [sol[1,end]-x1[1],sol[2,end]-x1[2],sol[3,end]-x1[3]]
end

V,k = 2.0, 2.0
x0,x1 = ones(3),zeros(3)
lambda0 = zeros(3)
Tf,N = 20.0, 50
dt = Tf/(N-1)

# initial guess 
X0,U0 = zeros(3,N), zeros(N-1)
for idx in 1:3
  X0[idx,:] = collect(linspace(x0[idx],x1[idx],N))
end

# SCP
constraints = Convex.Constraint[]
Xh,Uh = [X0], [U0]
max_iter,convergence_thresh = 50, 1e-3

for cvx_iter in 1:max_iter
  X,U = Variable(3,N), Variable(N-1)
  constraints += X[:,1] == x0
  constraints += X[:,end] == x1

  Jm = quadform(U,eye(N-1)) 
  for idx in 1:N-1
    f,Ak,Bk = dubins_dyn(X0[:,idx],U0[idx],dt) 
    constraints += X[:,idx+1] == X[:,idx] + dt*(f + Ak*(X[:,idx]-X0[:,idx]) + Bk*(U[idx]-U0[idx]))
  end

  prob = minimize(Jm,constraints)
  Convex.solve!(prob,SCSSolver())
  
  if prob.status != :Optimal 
    break
  end
 
  X0,U0 = copy(X.value), copy(vec(U.value))
  push!(Xh,X0)
  push!(Uh,U0)
  lambda0 = copy(constraints[1].dual)

  if norm(Xh[end]-Xh[end-1]) < convergence_thresh
    break
  end

  constraints = Convex.Constraint[]
  gc()
end

# PyPlot.figure()
# PyPlot.plot(Xh[end][1,:], Xh[end][2,:])
# PyPlot.quiver(X0[1,:], X0[2,:], V*cos.(X0[3,:]), V*sin.(X0[3,:]))

println("Guess for dual: $lambda0")
#nlsolve(shooting!, lambda0)
