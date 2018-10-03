%% Constants mass = 0.486;
mass = 6.0;
w_max = 0.5;
F_max = 2.5;
T_max = 2.5;

Tf = 25;
dt = 0.5;
N_mp = 120;






%% x = [qx qy qz qw wx wy wz] u = [Mx My Mz]
n = 6;  m = 3;
J = diag([0.3277, 0.3277, 0.5303]); % inertia tensor

f = @(x) [(x(1)*(2*x(1)*x(4) + 2*x(2)*x(5) + 2*x(3)*x(6)))/4 - (x(2)*x(6))/2 + (x(3)*x(5))/2 - (x(4)*(x(1)^2 + x(2)^2 + x(3)^2 - 1))/4;...
          (x(2)*(2*x(1)*x(4) + 2*x(2)*x(5) + 2*x(3)*x(6)))/4 + (x(1)*x(6))/2 - (x(3)*x(4))/2 - (x(5)*(x(1)^2 + x(2)^2 + x(3)^2 - 1))/4;...
          (x(3)*(2*x(1)*x(4) + 2*x(2)*x(5) + 2*x(3)*x(6)))/4 - (x(1)*x(5))/2 + (x(2)*x(4))/2 - (x(6)*(x(1)^2 + x(2)^2 + x(3)^2 - 1))/4;...
          (J(2,2)-J(3,3))*x(5)*x(6)/J(1,1);...
          (J(3,3)-J(1,1))*x(4)*x(6)/J(2,2);...
          (J(1,1)-J(2,2))*x(4)*x(5)/J(3,3)];

df = @(x) [(x(1)*x(4))/2+(x(2)*x(5))/2+(x(3)*x(6))/2, (x(1)*x(5))/2-x(6)/2-(x(2)*x(4))/2, x(5)/2+(x(1)*x(6))/2-(x(3)*x(4))/2, x(1)^2/4-x(2)^2/4-x(3)^2/4+1/4, x(3)/2+(x(1)*x(2))/2, (x(1)*x(3))/2-x(2)/2;...
        x(6)/2-(x(1)*x(5))/2+(x(2)*x(4))/2, (x(1)*x(4))/2+(x(2)*x(5))/2+(x(3)*x(6))/2, (x(2)*x(6))/2-x(4)/2-(x(3)*x(5))/2, (x(1)*x(2))/2-x(3)/2, -x(1)^2/4+x(2)^2/4-x(3)^2/4+1/4, x(1)/2+(x(2)*x(3))/2;...
        (x(3)*x(4))/2-(x(1)*x(6))/2-x(5)/2, x(4)/2-(x(2)*x(6))/2+(x(3)*x(5))/2, (x(1)*x(4))/2+(x(2)*x(5))/2+(x(3)*x(6))/2, x(2)/2+(x(1)*x(3))/2, (x(2)*x(3))/2-x(1)/2, -x(1)^2/4-x(2)^2/4+x(3)^2/4+1/4;...
        0,0,0, 0, (J(2,2)-J(3,3))*x(6)/J(1,1), (J(2,2)-J(3,3))*x(5)/J(1,1);...
        0,0,0, -(J(1,1)-J(3,3))*x(6)/J(2,2), 0, -(J(1,1)-J(3,3))*x(4)/J(2,2);...
        0,0,0, (J(1,1)-J(2,2))*x(5)/J(3,3), (J(1,1)-J(2,2))*x(4)/J(3,3), 0];
 
B = [zeros(3,3);...
    inv(J)];

% dummy function holders from q dynamics
xnorm = @(x) [0];
dxnorm = @(x) [0];

th_max = 359*pi/180;    th_tange = tan(th_max/4);
state_constr_low = [-th_tange*ones(3,1);-w_max/sqrt(3)*ones(3,1)];
state_constr = [state_constr_low, -state_constr_low];
ctrl_constr = T_max/sqrt(3)*[-ones(m,1) ones(m,1)];

x_eq = [quat2mrp([0;0;0;1]);0.01;-0.01;0.01];
u_eq = zeros(m,1);

% dynamics, quaternion norm, and terminal set constraint
c_L = [zeros((n)*(N_mp+1),1);0];
c_U = [zeros((n)*(N_mp+1),1);0.05];
c_nl = [c_L c_U];

test_state = [quat2mrp([sin(pi/4); 0; 0; cos(pi/4)]);-0.01;0.01;0];
filename = 'p_dynamics_mat.mat';





%% Dynamics and cost
Q = zeros(n); R = eye(m);

P = 15*eye(n);  %P(end,end) = 0;
alpha = 1e-3;
   
%% Misc parameters necessary
obs = struct('n_obs',0,'pos',[],'r',[]);
