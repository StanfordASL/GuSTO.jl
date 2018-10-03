function tau t2tau(t,tau0,tauf)
  % Converts from the interval [-1,1] to the interval [tau0,tauf]

  tau = zeros(size(t));
  for i = 1:numel(tau)
    tau(i) = 0.5*((tauf-tau0)*t(i) + (tauf + tau0));
  end
end
