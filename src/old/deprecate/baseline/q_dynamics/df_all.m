function out = df_all(df,x,n,N)
  out = zeros((N+1)*n);
  for j = 1:N+1
      out(1+(j-1)*n:j*n,1+(j-1)*n:j*n) = df(x(1+(j-1)*n:j*n));
  end
end
