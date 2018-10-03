function out = qnorm_all(qnorm,x,n,N)
  out = zeros((N+1),1);
  for j = 1:N+1
      out(j) = qnorm(x(1+n*(j-1):j*n));
  end
end
