function out = f_all(f,x,n,N)
  out = zeros((N+1)*n,1);
  for j = 1:N+1
      out(1+n*(j-1):j*n) = f(x(1+n*(j-1):j*n));
  end
end
