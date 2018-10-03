function phi = quat2angle(q)
  phi = 2*atan2(norm(q(1:4)), q(4));
end
