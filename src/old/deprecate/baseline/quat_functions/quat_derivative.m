function qdot = quat_derivative(q,w)
  omega_skew = [-skew(w), w; -w', 0];
  qdot = 0.5*omega_skew*q;
end
