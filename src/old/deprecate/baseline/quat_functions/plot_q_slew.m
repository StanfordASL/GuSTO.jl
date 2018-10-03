function [] = plot_q_slew(quat_matrix)
  % quat_matrix = 4xN matrix

  N = size(quat_matrix,2);
  X = eye(3);
  [Sx,Sy,Sz] = sphere(1e2);

  figure();
  for i = 1:N
    q = quat_matrix(:,i);
    dcm = quat2dcm(q);

    clf();
    hold on;
    for j = 1:3
      v = dcm*X(:,j);
      % surf(Sx,Sy,Sz);
      plot3([0, v(1)],...
            [0, v(2)],...
            [0, v(3)]);
      pause(0.1);
    end
  end
end
