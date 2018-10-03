function rotated_v = quat_rotate(q,v,vargin)
  % Eq. 183 in Shuster
  % v: 3 element vector to be rotated

  attitude_convention = "orientation";
  if nargin == 3
    attitude_convention = vargin{1}
  end

  qr = q;
  if attitude_convention == "rotation"
    qr = quat_inv(q);
  end

  v_bar = [v; 0];
  qprod = quat_multiply(qr, quat_multiply(v_bar, quat_inv(qr)));
  rotated_v = qprod(1:3);
end
