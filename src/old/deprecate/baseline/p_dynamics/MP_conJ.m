function conJ = MP_conJ(xu,Prob)%,D,n,m,N,df,B_full,Tf,P)
  n = Prob.user.n;
  N = Prob.user.N;
  m = Prob.user.m;

  conJ = zeros(n*(N+1)+1,(n+m)*(N+1));

  %% Dynamics constraints
  conJ(1:n*(N+1),1:n*(N+1)) = (2/Prob.user.Tf)*Prob.user.D - ...
              df_all(Prob.user.df,xu(1:n*(N+1)),n,N);

  conJ(1:n*(N+1),n*(N+1)+1:end) = -Prob.user.B_full;

  %% Terminal constraint
  xf = xu(n*N+1:n*(N+1));
  conJ(end,n*N+1:n*(N+1)) = Prob.user.P*(xf-Prob.user.x_eq);
end
