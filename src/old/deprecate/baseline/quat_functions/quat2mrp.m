function p = quat2mrp(q)
  p = q(1:3)/(1+q(4));
end
