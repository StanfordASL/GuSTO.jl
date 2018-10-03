function q_matrix = yershova_quaternions()
  % http://rotations.mitchell-lab.org/
  % Yershova, Anna, et al. "Generating uniform incremental grids on SO (3) using the Hopf fibration." 
  % The International journal of robotics research 29.7 (2010): 801-812.

  q_matrix = dlmread('../../quat_data.txt');
  for idx = 1:size(q_matrix,1)
    q_matrix(idx,:) = q_matrix(idx,:)./norm(q_matrix(idx,:));
  end
end
