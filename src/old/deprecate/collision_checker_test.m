close all; clc; clear all;

% Two Boxes
figure()
clearvars obstacles
obstacles{1} = [.1 .4 .4 .1 .1; .1 .1 .4 .4 .1];
obstacles{2} = [.4 .9 .9 .4 .4; .6 .6 .9 .9 .6];

[obstacles{1},K{1}] = rectangle_points([0.1;0.1;0.1], [0.4;0.4;0.4]);
[obstacles{2},K{2}] = rectangle_points([0.4;0.6;0.0], [0.9;0.9;0.9]);
test_pts = [0,0.25,0.4,0.6,0.6;...
            0,0.25,0.4,0.55,0.8;...
            0,0.25,0.4,0.4,0.4];

clf; hold on;
axis square; rectangle;
for i = 1:length(obstacles)
  % fill(obstacles{i}(1,:), obstacles{i}(2,:), 'r')
  trisurf(K{i}, obstacles{i}(1,:), obstacles{i}(2,:), obstacles{i}(3,:), 'FaceColor', 'red', 'Facealpha', 0.25);
end
plot3(test_pts(1,:), test_pts(2,:), test_pts(3,:), 'ob-', 'LineWidth',2, 'MarkerFaceColor','b');

% Windy Maze
figure();
clearvars obstacles
[obstacles{1}, K{1}] = rectangle_points([0;0.45;0],[0.9;0.55;0.01]);
[obstacles{2}, K{2}] = rectangle_points([0.25;0.2;0],[0.3;0.8;0.01]);
[obstacles{3}, K{3}] = rectangle_points([0.75;0.2;0],[0.8;0.8;0.01]);
[obstacles{4}, K{4}] = rectangle_points([0.5;0.7;0],[0.55;1;0.01]);
[obstacles{5}, K{5}] = rectangle_points([0.5;0;0],[0.55;0.3;0.01]);
            
test_pts = [.1,.25,.3,.5,.55,.75,.8,.9,.9,.8,.75,.55,.5,.3,.25,.15;...
          .65,.8,.8,.7,.7,.8,.8,.55,.45,.2,.2,.3,.3,.2,.2,.3;...
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
test_pts = [0.1, 0.275, 0.275, 0.275;...
            0.65, 0.8, 0.7, 0.7;...
            0, 0, 0, 0.01];

clf; hold on;
axis square; rectangle;
for i = 1:length(obstacles)
  % fill(obstacles{i}(1,:), obstacles{i}(2,:), 'r')
  trisurf(K{i}, obstacles{i}(1,:), obstacles{i}(2,:), obstacles{i}(3,:), 'FaceColor', 'red', 'Facealpha', 0.25);
end
plot3(test_pts(1,:), test_pts(2,:), test_pts(3,:), 'ob-', 'LineWidth',2, 'MarkerFaceColor','b');


function [pts,K] = rectangle_points(lo,hi)
  pts = [lo, hi];
  for i = 1:3
    v = lo;
    v(i) = hi(i);
    pts = [pts, v];
  end
  for i = 1:3
    v = hi;
    v(i) = lo(i);
    pts = [pts, v];
  end
  K = convhull(pts(1,:), pts(2,:), pts(3,:));
end
