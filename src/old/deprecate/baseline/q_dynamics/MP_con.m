function c = MP_con(xu,Prob)
  n = Prob.user.n;
  N = Prob.user.N;

  c = zeros((n+1)*(N+1)+1,1);

  %% Dynamics constraints
  c(1:n*(N+1)) = (2/Prob.user.Tf)*Prob.user.D*xu(1:n*(N+1)) -...
      f_all(Prob.user.f,xu(1:n*(N+1)),n,N);

  %% Quaternion norm constraints
  c(n*(N+1)+1:(n+1)*(N+1)) = qnorm_all(Prob.user.qnorm,xu(1:n*(N+1)),n,N);

  %% Terminal constraint
  xf = xu(n*N+1:n*(N+1));
  c(end) = (xf-Prob.user.x_eq)'*Prob.user.P*(xf-Prob.user.x_eq);
end
