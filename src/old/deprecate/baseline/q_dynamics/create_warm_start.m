load_PVTOL_config;
load MP_WARM_PVTOL.mat;

ts = 0.5*mp_warm.Tp*(mp_warm.s_t+1);
N = numel(mp_warm.s_t)-1;

%% Double integrator dynamics
A = [zeros(3) eye(3); zeros(3,6)];
B = [zeros(3); eye(3)];
C = eye(6);
D = zeros(6,3);
system = ss(A,B,C,D);

ts = linspace(0,mp_warm.Tp,N+1);
ctrl = ones(m,N+1);
output = lsim(system, mp_warm.ctrl', ts, test_state);
mp_warm.state = output;

for i = 1:(N+1)
    mp_warm.result.x_0((n+m)*(i-1)+1:i*(n+m)) = [output(i,:)'; ctrl(:,i)];
    mp_warm.result.x_k((n+m)*(i-1)+1:i*(n+m)) = [output(i,:)'; ctrl(:,i)];
    mp_warm.result.SOL.xs((n+m)*(i-1)+1:i*(n+m)) = [output(i,:)'; ctrl(:,i)];
end

save('MP_WARM_PVTOL.mat','mp_warm');

%% Quaternion dynamics
% x0 = [-0.5  0.5 0.5 -0.5 0.1 0.1 0.1]';
% J = diag([10 5 7]);
% Ts = repmat([3 1 -7]', 1,numel(ts));
% Xs = zeros(7,numel(ts));
% Xs(:,1) = x0;
% for i = 2:numel(ts)
%    out =  ode45(@q_dynamics,[ts(i-1) ts(i)],Xs(:,i-1),[],Ts(:,i),J);
%    q_out = out.y(1:4,end);
%    w_out = out.y(5:7,end);
%    Xs(:,i) = [q_out./norm(q_out); w_out];
% end
% 
% mp_warm.solve = 1;
% mp_warm.state = Xs';
% save('MP_WARM_PVTOL.mat', 'mp_warm');

%% MRP dynamics

%% Functions

function xdot = q_dynamics(t,x,T,J)
    % xdot = dynamics(x,T,J)
    %
    % INPUTS:
    %   x = 7x1 vector [q,w]
    %   T = 3x1 torque command
    %   J = 3x3 robot inertia tensor
    %
    % Eq. 305 in Shuster
    % q_dot = 0.5 * omega_skew * q;

    q = x(1:4,:);
    w = x(5:7,:);

    q_dot = quat_derivative(q,w);
    w_dot = inv(J)*(T - cross(w, J*w));

    xdot = [q_dot; w_dot];
end

function f = p_dynamics(x,T,J)
    % f = dynamics(x,T,J)
    %
    % INPUTS:
    %   x = 6x1 vector [p,w]
    %   T = 3x1 torque command
    %   J = 3x3 robot inertia tensor

    p = x(1:3,:);
    w = x(4:6,:);

    p_dot = mrp_derivative(p,w);
    w_dot = inv(J)*(T - cross(w, J*w));

    f = [p_dot; w_dot];
end
