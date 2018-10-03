function qinv = quat_inv(q)
  % Eq. 177 in Shuster
  qinv = q; 
  qinv(1:3) = -1*qinv(1:3);
  qinv = qinv./norm(qinv);
end
