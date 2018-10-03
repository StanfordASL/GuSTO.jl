function rb = robot()
  rb = struct();
  rb.mass = 6.0;
  rb.n_thrusters = 12;

  s = 0.5*0.305; % each side of cube is 30.5cm
  rb.s = s;
  rb.r = sqrt(3)*s;   % inflate to sphere

  rb.J = diag([0.3277, 0.3277, 0.5303]);
  rb.hard_limit_vel   = 0.0200;
  rb.hard_limit_accel = 0.0015;
  rb.hard_limit_omega = 0.0750;
  rb.hard_limit_alpha = 0.0250;
  
  % thrusters:
  % 1. RX+    2. LX+      3. RX-      4. LX-
  % 5. AY+    6. FY+      7. AY-      8. FY-
  % 9. LZ+    10. RZ+     11. LZ-     12. RZ-

  % guesses based upon pictures in paper
  rb.thruster_normals = [1., 0, 0;, 1, 0, 0; -1, 0, 0; -1, 0, 0; 0, 1, 0; 0, 1, 0; 0, -1, 0; 0, -1, 0; 0, 0, 1; 0, 0, 1; 0, 0, -1; 0, 0, -1];
  rb.thruster_positions = [s, s, -s; s, -s, s; -s, s, s; -s, -s, -s; -s, s, s; s, s, -s; -s, -s, -s; s, -s, s; s, -s, s; -s, s, s; -s, -s, -s; s, s, -s];

  % thruster allocation matrix
  rb.G = zeros(6, rb.n_thrusters);
  for idx = 1:rb.n_thrusters
    rb.G(1:3,idx) = rb.thruster_normals(idx,:);
    rb.G(4:6,idx) = cross(rb.thruster_positions(idx,:), rb.thruster_normals(idx,:));
  end
end
