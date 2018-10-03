function taus = t2tau(ts,tau0,tauf)
  % Converts from the interval [-1,1] to the interval [tau0,tauf]

  taus = zeros(size(ts));
  for i = 1:numel(tau)
    t = ts(i);
    taus(t) = 0.5*((1-t)*tau0 + (1+t)*tauf);
  end
end
