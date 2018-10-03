function ts = tau2t(taus)
  % Converts from the interval [tau0,tauf] to the interval [-1,1]

  ts = zeros(size(taus));
  tau0 = taus(1);
  tauf = taus(end);
  Delta = tauf-tau0;

  for i = 1:numel(taus)
    ts(i) = (2*taus(i)-Delta)/Delta;
  end
end
