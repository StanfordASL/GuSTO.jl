function dcm = quat2dcm(q, vargin)
  % Eq. 97 in Shuster (note sin term sign flipped based on skew matrix convention)
  % Eq. 1.26 on p20 in Spacecraft Dynamics and Control by Anton de Ruiter

  attitude_convention = "rotation";
  if nargin == 2
    attitude_convention = vargin{1};
  end

  phi = quat2angle(q);
  n = quat2vec(q);
  dcm = cos(phi)*eye(3) + (1-cos(phi))*n*n' - sin(phi)*skew(n);

  if attitude_convention == "rotation"
    dcm = dcm';
  elseif attitude_convention ~= "orientation"
    warn("Must specify if attitude is rotation or orientation! Defaulting to orientation.")
  end
end
