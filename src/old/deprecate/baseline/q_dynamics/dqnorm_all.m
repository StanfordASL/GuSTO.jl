function out = dqnorm_all(qnormJ,x,n,N)
  out = zeros(N+1,n*(N+1));
  for j = 1:N+1
      out(j,1+n*(j-1):j*n) = qnormJ(x(1+n*(j-1):j*n));
  end
end
