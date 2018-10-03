function a = quat2vec(q)
  % canonicalize: true corresponds to positive scalar component q for dual representation
  if norm(q(1:3)) == 0
    a = zeros(3,1);
    return
  end

  a = sign(q(4))*q(1:3)/norm(q(1:3));
end
