%% x = [qx qy qz qw wx wy wz] u = [Mx My Mz]  u=Gf 

clear all; close all; clc;
warning('off','MATLAB:nargchk:deprecated');
addpath('../quat_functions');

mass = 6.0;
w_max = 0.5;
F_max = 2.5;
T_max = 2.5;

Tf = 25;
dt = 0.01;
N = 120; %Number of collocation points

% load('scp_initial_guess.mat'); % SBMP intermediate waypoints 
filename = 'q_dynamics_mat.mat';


%% Dynamics
n = 7;  m = 3;  m_f = 12; 
J = 0.1083*eye(3); % inertia tensor

f = @(x)[(x(2)*x(7))/2 - (x(3)*x(6))/2 + (x(4)*x(5))/2;...
   (x(3)*x(5))/2 - (x(1)*x(7))/2 + (x(4)*x(6))/2;...
   (x(1)*x(6))/2 - (x(2)*x(5))/2 + (x(4)*x(7))/2;...
 - (x(1)*x(5))/2 - (x(2)*x(6))/2 - (x(3)*x(7))/2;...
    (J(2,2)*x(6)*x(7) - J(3,3)*x(6)*x(7))/J(1,1);...
    (-J(1,1)*x(5)*x(7) + J(3,3)*x(5)*x(7))/J(2,2);...
    (J(1,1)*x(5)*x(6) - J(2,2)*x(5)*x(6))/J(3,3)];

df = @(x) [0,  x(7)/2, -x(6)/2, x(5)/2,  x(4)/2, -x(3)/2,  x(2)/2;...
  -x(7)/2,     0,  x(5)/2, x(6)/2,  x(3)/2,  x(4)/2, -x(1)/2;...
   x(6)/2, -x(5)/2,     0, x(7)/2, -x(2)/2,  x(1)/2,  x(4)/2;...
  -x(5)/2, -x(6)/2, -x(7)/2,    0, -x(1)/2, -x(2)/2, -x(3)/2;...
  0, 0, 0, 0, 0, (J(2,2)*x(7) - J(3,3)*x(7))/J(1,1),  (J(2,2)*x(6) - J(3,3)*x(6))/J(1,1);...
  0, 0, 0, 0, -(J(1,1)*x(7) - J(3,3)*x(7))/J(2,2), 0, -(J(1,1)*x(5) - J(3,3)*x(5))/J(2,2);...
  0, 0, 0, 0,  (J(1,1)*x(6) - J(2,2)*x(6))/J(3,3), (J(1,1)*x(5) - J(2,2)*x(5))/J(3,3), 0];

B = [zeros(4,3);...
    inv(J)];

qnorm = @(x) [sqrt(x(1)^2+x(2)^2+x(3)^2+x(4)^2)]-1;
dqnorm = @(x) [x(1)/sqrt(x(1)^2+x(2)^2+x(3)^2+x(4)^2);...
              x(2)/sqrt(x(1)^2+x(2)^2+x(3)^2+x(4)^2);...
              x(3)/sqrt(x(1)^2+x(2)^2+x(3)^2+x(4)^2);...
              x(4)/sqrt(x(1)^2+x(2)^2+x(3)^2+x(4)^2);...
              zeros(3,1)]';


% Thruster allocation matrix
s = 0.5*0.305;
n_thrusters = 12;
thruster_normals    = [1,0,0; 1,0,0; -1,0,0; -1,0,0; 0,1,0; 0,1,0; 0,-1,0; 0,-1,0; 0,0,1; 0,0,1; 0,0,-1; 0,0,-1];
thruster_positions  = [s,s,-s; s,-s,s; -s,s,s; -s,-s,-s; -s,s,s; s,s,-s; -s,-s,-s; s,-s,s; s,-s,s; -s,s,s; -s,-s,-s; s,s,-s];
G = zeros(3, n_thrusters);
for idx  = 1:n_thrusters
 % G(1:3,idx) = thruster_normals(idx,:);
 G(1:3,idx) = cross(thruster_positions(idx,:), thruster_normals(idx,:)); 
end

%% Constraints
% State bounds
x_L = -[ones(4,1);w_max/sqrt(3)*ones(3,1)];
x_U = [ones(4,1);w_max/sqrt(3)*ones(3,1)];
state_constr = [x_L x_U];

% Control bounds
u_L = -T_max/sqrt(3)*ones(m,1);
u_U = T_max/sqrt(3)*ones(m,1);
ctrl_constr = [u_L u_U];
f_L = zeros(m_f,1);
f_U = 1e10*ones(m_f,1);
           
x_eq = [0;0;0;-1;0.01;-0.01;0.01];
u_eq = zeros(m,1);
f_eq = zeros(m_f,1);

% linear equality constraints
A = zeros(m*(N+1),(n+m+m_f)*(N+1));
A(:,n*(N+1)+1:(n+m)*(N+1)) = kron(eye(N+1),eye(m)); 
A(:,(n+m)*(N+1)+1:(n+m+m_f)*(N+1)) = kron(eye(N+1),-G);
b_L = zeros(m*(N+1),1);
b_U = zeros(m*(N+1),1);

% nonlinear equality constraints
% dynamics, quaternion norm, and terminal set constraint
c_L = [zeros((n+1)*(N+1),1);0];
c_U = [zeros((n+1)*(N+1),1);0.05];

test_state = [sin(pi/4); 0; 0; cos(pi/4);-0.01;0.01;0];


% Dynamics and cost
Q = zeros(n); R = eye(m+m_f);

P = 15*eye(n);  %P(end,end) = 0;
alpha = 1e-3;

%% Misc parameters necessary
obs = struct('n_obs',0,'pos',[],'r',[]);
  


%% Setup problem
% CGL nodes
[s_t,w] = clencurt(N); %t_t: [-1, 1] : <-> : [0, Tf]
s = fliplr(s_t); %t: [1, -1]

% Final solution interpolation matrix
tau_full = 0:dt:Tf; 
s_e_full = (2*tau_full - Tf)/Tf; %[-1, 1]

% idxs = zeros(numel(ts_wps),1);
% for k = 1:numel(ts_wps)
%   t_wp = ts_wps(k);
%   [val,idxs(k)] = min(abs(s_t-t_wp)); % TODO(acauligi): check if s_t or s here
% end

% Lagrange polynomial evaluation at the interpolation points
L_e = compute_Lagrange(length(s_e_full)-1,N,s_e_full,s_t);

% Get Differentiation matrix
D = ChebyshevDiffMatrix(N,s); %arranged for forward time
D = kron(D,eye(n));

x_eq_all = kron(ones(N+1,1),x_eq);
u_eq_all = kron(ones(N+1,1),u_eq);
f_eq_all = kron(ones(N+1,1),f_eq);

xu_eq = [x_eq_all;u_eq_all;f_eq_all];

Q_bar = kron(diag(w),Q); R_bar = kron(diag(w),R);
Q_tilde = Q_bar + kron(diag([zeros(N,1);(2/Tf)]),P);

F = blkdiag(Q_tilde,R_bar);

B_full = kron(eye(N+1),B);

xu_L = [kron(ones(N+1,1),x_L);
        kron(ones(N+1,1),u_L)
        kron(ones(N+1,1),f_L)];
xu_U = [kron(ones(N+1,1),x_U);
        kron(ones(N+1,1),u_U); 
        kron(ones(N+1,1),f_U)];
    
% for k = 1:numel(idxs)
%   idx = idxs(k);
%   xu_L(n*(idx-1)+1:n*(idx-1)+4) = wps(1:4,k);
%   xu_U(n*(idx-1)+1:n*(idx-1)+4) = wps(1:4,k);
% end
       
MPC_cost = @(xu) (Tf/2)*(xu-xu_eq)'*F*(xu-xu_eq);
MPC_grad = @(xu) Tf*F*(xu-xu_eq);
MPC_hess = @(xu) Tf*F;

xu0 = zeros((n+m+m_f)*(N+1),1);

Prob = conAssign(MPC_cost,MPC_grad,MPC_hess,[],...
            xu_L,xu_U,'MP', xu0,...
            [], 0, A, b_U, b_L,...
            'MP_con','MP_conJ',[],[],...
            c_L,c_U,...
            [],[],[],[]);

Prob.user.x_act     = zeros(n,1);
Prob.user.D         = D;
Prob.user.n         = n;
Prob.user.m         = m;
Prob.user.N         = N;
Prob.user.f         = f;
Prob.user.df        = df;
Prob.user.B_full    = B_full;
Prob.user.qnorm     = qnorm;
Prob.user.dqnorm    = dqnorm;
Prob.user.Tf        = Tf;
Prob.user.P         = P;
Prob.user.x_eq      = x_eq;

%% Load warm start
if exist(filename, 'file') == 2
    load(filename);
else
    mp_warm = struct('Tp',Tf,'shift',0,'sol',0,...
                's_t',s_t,'state',[],'ctrl',[],'result',[]);
end


%% Solve problem
% adjust initial nominal state guess
for i = 1:n
    if test_state(i) > state_constr(i,2)
        test_state(i) = state_constr(i,2)-0.01*(state_constr(i,2)-state_const(i,1));
    elseif test_state(i) < state_constr(i,1)
        test_state(i) = state_constr(i,1)+0.01*(state_constr(i,2)-state_const(i,1));
    end
end

% Solution guess
if (mp_warm.sol==0) %got nothing
    In = eye(n);
    x0 = zeros((N+1)*n,1);
    for i = 1:n
        x0 = x0 + kron(linspace(test_state(i),x_eq(i),N+1), In(i,:))';
    end
    
    taus = cos(pi*(0:N)/N);
    q1 = [test_state(4), test_state(1:3)'];
    q2 = [x_eq(4), x_eq(1:3)'];
    for idx = 1:N+1
      t_slew = (taus(idx)+1)/2;
      qout = quatinterp(q1,q2,t_slew,'slerp');
      x0(n*(idx-1)+1:n*(idx-1)+4) = [qout(2);qout(3);qout(4);qout(1)];
    end
    
    Im = eye(m);
    u0 = zeros((N+1)*m,1);
    for j = 1:m
        u0 = u0 + kron(u_eq(j)*ones(N+1,1), Im(:,j));
    end

    Im = eye(m_f);
    f0 = zeros((N+1)*m_f,1);
    for j = 1:m
        f0 = f0 + kron(f_eq(j)*ones(N+1,1), Im(:,j));
    end
    
    Prob = modify_x_0(Prob,[x0;u0;f0]);
elseif (mp_warm.sol==0.5) % have some guess of homotopy
    Prob = modify_x_0(Prob,mp_warm.result.x_k);
end

% Update constraint information
Prob.user.x_act = test_state; 

% Recall warm solution
if (mp_warm.sol==1) %have actual full solution
    Prob = WarmDefSOL('snopt',Prob,mp_warm.result);
end

Prob = ProbCheck(Prob,'snopt');

%% Solve
tic();
Result = snoptTL(Prob);
toc()

converged = Result.Inform; %GOOD: {1,2,3}

% Compute trajectories
MP_state = zeros(size(L_e,2),n);
x_nom = zeros(N+1,n);
for i = 1:n
    c = Result.x_k(i:n:n*(N+1)-(n-i))';
    MP_state(:,i) = (c*L_e)';
    x_nom(:,i) = c';
end


MP_ctrl = zeros(size(L_e,2),m);
u_nom = zeros(N+1,m);
for j = 1:m
    c = Result.x_k(n*(N+1)+j:m:end-(m-j))';
    MP_ctrl(:,j) = (c*L_e)';
    u_nom(:,j) = c';
end

disp('MP:'); disp(converged);
