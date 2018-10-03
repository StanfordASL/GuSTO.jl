function qprod = quat_multiply(q1,q2)
  % Eq. 169-171 in Shuster
  % q1 * q2 <==> R(q1) * R(q2)
  qprod = zeros(4,1);
  qprod(1:3) = q1(4)*q2(1:3) + q2(4)*q1(1:3) - cross(q1(1:3), q2(1:3));
  qprod(4) = q1(4)*q2(4) - dot(q1(1:3), q2(1:3));
end
