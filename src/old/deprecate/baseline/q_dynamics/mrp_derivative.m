function pdot = mrp_derivative(p,w)
  pdot = 0.25*((1-norm(p)^2)*eye(3) +2*skew(p) + 2*p*p')*w;
end
