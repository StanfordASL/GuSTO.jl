function dcm = mrp2dcm(p)
  % Eq. 255b in Shuster
  % NOTE (6/26): switched +4 term to -4 term because that's what I 
  % found everywhere else and unit tests pass with this only

  dcm = eye(3) + (8*skew(p)^2 - 4*(1 - norm(p)^2)*skew(p))/(1+norm(p)^2)^2;

end


function skew_v = skew(v)
    skew_v = [0 -v(3) v(2) ; v(3) 0 -v(1) ; -v(2) v(1) 0 ];
end