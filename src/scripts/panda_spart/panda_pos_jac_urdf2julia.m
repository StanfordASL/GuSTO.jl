%% Computes the position and position jacobian of each joint
%   given a urdf description of a robot
%   1st edit: January 19, 2019
%% --- Clean and clear ---%
clc
close all
clear

%% --- URDF filename ---%
filename='panda.urdf';

%% --- Create robot model ---%
[robot,robot_keys] = urdf2robot(filename);

%% --- Parameters ---%
R0=eye(3);
r0=zeros(3,1);
syms q1 q2 q3 q4 q5 q6 q7
assume(q1, 'real');
assume(q2, 'real');
assume(q3, 'real');
assume(q4, 'real');
assume(q5, 'real');
assume(q6, 'real');
assume(q7, 'real');
qm = [q1;q2;q3;q4;q5;q6;q7];

%--- Kinematics ---%
%% Kinematics
[RJ,RL,rJ,rL,e,g]=Kinematics(R0,r0,qm,robot);
r_EE = simplify(rJ(:,end), 'Steps', 200);
J_p_EE = jacobian(rJ(:,end), qm);

Q1 = zeros(robot.n_q,1);
Q2 = [pi/4; pi/2; 0; -0.8*pi; 0.25; 0.75*pi; 0;];
Q3 = [-pi/4; pi/2; 0; -0.8*pi; 0.25; 0.75*pi; 0];
Q4 = [pi/4; 0; 0; -0.8*pi; 0.25; 0.75*pi; 0];

eval(subs(r_EE, qm, Q1));
eval(subs(J_p_EE, qm, Q1));

eval(subs(r_EE, qm, Q2));
eval(subs(r_EE, qm, Q3));

%% Symbolic part - computation of positions
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16
assume(x1, 'real');
assume(x2, 'real');
assume(x3, 'real');
assume(x4, 'real');
assume(x5, 'real');
assume(x6, 'real');
assume(x7, 'real');
assume(x8, 'real');
assume(x9, 'real');
assume(x10, 'real');
assume(x11, 'real');
assume(x12, 'real');
assume(x13, 'real');
assume(x14, 'real');
embedding = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14];

% % Replace half cosines and sines
% % sin(a/2)^2 = (1-cos(a))/2
% % cos(a/2)^2 = (1+cos(a))/2
% % sin(a/2) cos(a/2) = sin(a)/2
% 
% % p_EE = rJ(:,end);
% rJ = subs(rJ, sin(q1/2)*cos(q1/2), sin(q1)/2);
% rJ = subs(rJ, cos(q1/2)*sin(q1/2), sin(q1)/2);
% rJ = subs(rJ, sin(q1/2)^2, (1-cos(q1))/2);
% rJ = subs(rJ, cos(q1/2)^2, (1+cos(q1))/2);
% rJ = subs(rJ, cos(q1), x1);
% rJ = subs(rJ, sin(q1), x2);
% 
% rJ = subs(rJ, sin(q2/2)*cos(q2/2), sin(q2)/2);
% rJ = subs(rJ, cos(q2/2)*sin(q2/2), sin(q2)/2);
% rJ = subs(rJ, sin(q2/2)^2, (1-cos(q2))/2);
% rJ = subs(rJ, cos(q2/2)^2, (1+cos(q2))/2);
% rJ = subs(rJ, cos(q2), x3);
% rJ = subs(rJ, sin(q2), x4);
% 
% rJ = subs(rJ, sin(q3/2)*cos(q3/2), sin(q3)/2);
% rJ = subs(rJ, cos(q3/2)*sin(q3/2), sin(q3)/2);
% rJ = subs(rJ, sin(q3/2)^2, (1-cos(q3))/2);
% rJ = subs(rJ, cos(q3/2)^2, (1+cos(q3))/2);
% rJ = subs(rJ, cos(q3), x5);
% rJ = subs(rJ, sin(q3), x6);
% 
% rJ = subs(rJ, sin(q4/2)*cos(q4/2), sin(q4)/2);
% rJ = subs(rJ, cos(q4/2)*sin(q4/2), sin(q4)/2);
% rJ = subs(rJ, sin(q4/2)^2, (1-cos(q4))/2);
% rJ = subs(rJ, cos(q4/2)^2, (1+cos(q4))/2);
% rJ = subs(rJ, cos(q4), x7);
% rJ = subs(rJ, sin(q4), x8);
% 
% rJ = subs(rJ, sin(q5/2)*cos(q5/2), sin(q5)/2);
% rJ = subs(rJ, cos(q5/2)*sin(q5/2), sin(q5)/2);
% rJ = subs(rJ, sin(q5/2)^2, (1-cos(q5))/2);
% rJ = subs(rJ, cos(q5/2)^2, (1+cos(q5))/2);
% rJ = subs(rJ, cos(q5), x9);
% rJ = subs(rJ, sin(q5), x10);
% 
% rJ = subs(rJ, sin(q6/2)*cos(q6/2), sin(q6)/2);
% rJ = subs(rJ, cos(q6/2)*sin(q6/2), sin(q6)/2);
% rJ = subs(rJ, sin(q6/2)^2, (1-cos(q6))/2);
% rJ = subs(rJ, cos(q6/2)^2, (1+cos(q6))/2);
% rJ = subs(rJ, cos(q6), x11);
% rJ = subs(rJ, sin(q6), x12);
% 
% rJ = subs(rJ, sin(q7/2)*cos(q7/2), sin(q7)/2);
% rJ = subs(rJ, cos(q7/2)*sin(q7/2), sin(q7)/2);
% rJ = subs(rJ, sin(q7/2)^2, (1-cos(q7))/2);
% rJ = subs(rJ, cos(q7/2)^2, (1+cos(q7))/2);
% rJ = subs(rJ, cos(q7), x13);
% rJ = subs(rJ, sin(q7), x14);
% 
% rJ = simplify(rJ, 'Steps', 200);
% 
% disp('Finished computing all positions!')
% ccode(rJ,   'file','p_EE_panda');

%% Compute Jacobians and export everything to c-code
% for i = 1 : robot.n_q
%     fprintf('Computing and exporting position and Jac. of joint #%d\n', i)
%     
%     p_joint_i = rJ(:, i);
%     Jp_joint_i = simplify(jacobian(p_joint_i, embedding), 'Steps', 200);
%     
%     filename = sprintf('p_panda_joint_%d.jl',  i);
%     ccode(p_joint_i,  'file', filename);
%     position_or_jacobian_C_to_julia(filename);
%     
%     filename = sprintf('Jp_panda_joint_%d.jl',  i);
%     ccode(Jp_joint_i, 'file', filename);
%     position_or_jacobian_C_to_julia(filename);
% end


z_hat_I = [0;0;1];
x_hat_EE = [1;0;0];
R_I2EE = RJ(:,:,end);

pointing_constraint   = simplify(1 - dot(x_hat_EE, R_I2EE*z_hat_I), 'Steps', 200);
pointing_constraint = subs(pointing_constraint, sin(q1/2)*cos(q1/2), sin(q1)/2);
pointing_constraint = subs(pointing_constraint, cos(q1/2)*sin(q1/2), sin(q1)/2);
pointing_constraint = subs(pointing_constraint, sin(q1/2)^2, (1-cos(q1))/2);
pointing_constraint = subs(pointing_constraint, cos(q1/2)^2, (1+cos(q1))/2);
pointing_constraint = subs(pointing_constraint, cos(q1), x1);
pointing_constraint = subs(pointing_constraint, sin(q1), x2);

pointing_constraint = subs(pointing_constraint, sin(q2/2)*cos(q2/2), sin(q2)/2);
pointing_constraint = subs(pointing_constraint, cos(q2/2)*sin(q2/2), sin(q2)/2);
pointing_constraint = subs(pointing_constraint, sin(q2/2)^2, (1-cos(q2))/2);
pointing_constraint = subs(pointing_constraint, cos(q2/2)^2, (1+cos(q2))/2);
pointing_constraint = subs(pointing_constraint, cos(q2), x3);
pointing_constraint = subs(pointing_constraint, sin(q2), x4);

pointing_constraint = subs(pointing_constraint, sin(q3/2)*cos(q3/2), sin(q3)/2);
pointing_constraint = subs(pointing_constraint, cos(q3/2)*sin(q3/2), sin(q3)/2);
pointing_constraint = subs(pointing_constraint, sin(q3/2)^2, (1-cos(q3))/2);
pointing_constraint = subs(pointing_constraint, cos(q3/2)^2, (1+cos(q3))/2);
pointing_constraint = subs(pointing_constraint, cos(q3), x5);
pointing_constraint = subs(pointing_constraint, sin(q3), x6);

pointing_constraint = subs(pointing_constraint, sin(q4/2)*cos(q4/2), sin(q4)/2);
pointing_constraint = subs(pointing_constraint, cos(q4/2)*sin(q4/2), sin(q4)/2);
pointing_constraint = subs(pointing_constraint, sin(q4/2)^2, (1-cos(q4))/2);
pointing_constraint = subs(pointing_constraint, cos(q4/2)^2, (1+cos(q4))/2);
pointing_constraint = subs(pointing_constraint, cos(q4), x7);
pointing_constraint = subs(pointing_constraint, sin(q4), x8);

pointing_constraint = subs(pointing_constraint, sin(q5/2)*cos(q5/2), sin(q5)/2);
pointing_constraint = subs(pointing_constraint, cos(q5/2)*sin(q5/2), sin(q5)/2);
pointing_constraint = subs(pointing_constraint, sin(q5/2)^2, (1-cos(q5))/2);
pointing_constraint = subs(pointing_constraint, cos(q5/2)^2, (1+cos(q5))/2);
pointing_constraint = subs(pointing_constraint, cos(q5), x9);
pointing_constraint = subs(pointing_constraint, sin(q5), x10);

pointing_constraint = subs(pointing_constraint, sin(q6/2)*cos(q6/2), sin(q6)/2);
pointing_constraint = subs(pointing_constraint, cos(q6/2)*sin(q6/2), sin(q6)/2);
pointing_constraint = subs(pointing_constraint, sin(q6/2)^2, (1-cos(q6))/2);
pointing_constraint = subs(pointing_constraint, cos(q6/2)^2, (1+cos(q6))/2);
pointing_constraint = subs(pointing_constraint, cos(q6), x11);
pointing_constraint = subs(pointing_constraint, sin(q6), x12);

pointing_constraint = subs(pointing_constraint, sin(q7/2)*cos(q7/2), sin(q7)/2);
pointing_constraint = subs(pointing_constraint, cos(q7/2)*sin(q7/2), sin(q7)/2);
pointing_constraint = subs(pointing_constraint, sin(q7/2)^2, (1-cos(q7))/2);
pointing_constraint = subs(pointing_constraint, cos(q7/2)^2, (1+cos(q7))/2);
pointing_constraint = subs(pointing_constraint, cos(q7), x13);
pointing_constraint = subs(pointing_constraint, sin(q7), x14);

pointing_constraint = simplify(pointing_constraint, 'Steps', 200);

J_pointing_constraint = simplify(jacobian(pointing_constraint, embedding), 'Steps', 200);
H_pointing_constraint = simplify(jacobian(J_pointing_constraint, embedding), 'Steps', 200);

filename = 'pointing_constraint.jl';
ccode(pointing_constraint,  'file', filename);
position_or_jacobian_C_to_julia(filename);

filename = 'J_pointing_constraint.jl';
ccode(J_pointing_constraint,  'file', filename);
position_or_jacobian_C_to_julia(filename);

filename = 'H_pointing_constraint.jl';
ccode(H_pointing_constraint, 'file', filename);
position_or_jacobian_C_to_julia(filename);
