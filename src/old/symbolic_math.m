%% Discretization of double integrator dynamics
% Fundamentals of Linear State Space Systems (John Bay)
syms t t0 tau
assume(t, 'real');
assume(t0, 'real');
assume(tau, 'real');
X0 = sym('X0', [6,1], 'real');
U = sym('U', [3,1], 'real');
A = [zeros(3,3), eye(3); zeros(3,6)];
B = [zeros(3,3); eye(3)];

X = expm(A*(t-t0))*X0 + int(expm(A*(t-tau))*B*U, tau, 0, t);    % Eq. 6.6






%% Kinodynamic RRT*: Asymptotically optimal motion planning for robots with linear dynamics - Webb & van den Berg, 2013
% double integrator dynamics and c = 0
clear all; clc;
syms t tau tp r1 r2 r3 m
assume(t, 'real');
assume(tau, 'real');
assume(tp, 'real');
assume(r1, 'real');
assume(r2, 'real');
assume(r3, 'real');
assume(m, 'real');
A = [0,1; 0,0];     A = kron(A, eye(3));
B = [0;1];          B = kron(B, eye(3));
R = diag([r1,r2,r3]);

% Eq. 6
G_integrand = expm(A*(t - tp))*B*inv(R)*B'*expm(A'*(t-tp));
G = int(G_integrand,tp,0,t);
G = subs(G, t, tau);

% Eq. 11
xi = sym('xi', [6,1], 'real');
xf = sym('xf', [6,1], 'real');

xbar = expm(A*tau)*xi;  % Eq. 8

c = tau + (xf - xbar)'*inv(G)*(xf - xbar);

% Eq. 14
d = inv(G)*(xf - xbar);

% Eq. 11
cdot = 1 - 2*(A*xf)'*d - d'*B*inv(R)*B'*d;

% Eq. 20
soln = expm([A, B*inv(R)*B'; zeros(size(A)), -A'] * (t - tau))*[xf; d];
xt = soln(1:6);
yt = soln(7:12);
u = inv(R)*B'*yt;






%% A Survey of Attitude Representations - Malcom D. Shuster (1993)
clear all; clc;
syms Jxx Jyy Jzz Jxy Jxz Jyz
assume(Jxx, 'real');
assume(Jyy, 'real');
assume(Jzz, 'real');
p = sym('p', [3,1], 'real');
wi = sym('wi', [3,1], 'real');
wb = sym('wb', [3,1], 'real');
M = sym('M', [3,1], 'real');
J = diag([Jxx,Jyy,Jzz]);

pb_dot = 0.25*((1-p'*p)*wb + 2*cross(p,wb) + 2*dot(wb,p)*p);  % Eq. 338
wbdot = inv(J)*(M - cross(wb,J*wb));
Apb = [jacobian(pb_dot, [p;wb]); ...
    jacobian(wbdot, [p;wb])];
Bpb = [jacobian(pb_dot, [M]); ...
    jacobian(wbdot, [M])];

pi_dot = 0.25*((1-p'*p)*wi - 2*cross(p,wi) + 2*dot(wi,p)*p);  % Eq. 340
widot = inv(J)*(M - cross(wi,J*wi));
Api = [jacobian(pi_dot, [p;wi]); ...
    jacobian(widot, [p;wi])];
Bpi = [jacobian(pi_dot, [M]); ...
    jacobian(widot, [M])];





%% Rigid body dynamics
clear all; clc;
q = sym('q', [4,1], 'real');
w = sym('w', [3,1], 'real');
M = sym('M', [3,1], 'real');
%J = sym('J', [3,3], 'real');
syms Jxx Jyy Jzz Jxy Jxz Jyz
assume(Jxx, 'real');
assume(Jyy, 'real');
assume(Jzz, 'real');
%assume(Jxy, 'real');
%assume(Jyz, 'real');
%assume(Jxz, 'real');
%J = [Jxx,Jxy,Jxz;...
%    Jxy,Jyy,Jyz;...
%    Jxz,Jyz,Jzz];
J = diag([Jxx,Jyy,Jzz]);

omega_skew = [0, w(3),-w(2),w(1);...
            -w(3),0,w(1),w(2);...
            w(2),-w(1),0,w(3);...
            -w(1),-w(2),-w(3),0];
qdot = 0.5 * omega_skew * q;
wdot = inv(J)*(M - cross(w,J*w));

x = [q;w];
xdot = [qdot;wdot];
A = jacobian(xdot, x);
B = jacobian(xdot, M);
Hess = cell(numel(x),1);
for k = 1:numel(x)
    Hess{k} = simplify(hessian(xdot(k),[x;M]));
end




%% Motion planning with sequential convex optimization and convex collision checking- Schulman et al. (2014)
% evaluating Eq. 24 for "stop-and-turn" needle strategy
clear all; clc;
syms phi D Dk
assume(phi, 'real');
assume(D, 'real');
assume(Dk, 'real');

% exp map on SE(3) from A Mathematical Introduction to Robotic Manipulation
v = sym('v', [3,1], 'real'); % coefficients of translational se(3) matrices
w = sym('w', [3,1], 'real'); % coefficients of rotational se(3) matrices
se3 = [v;w];
th = sqrt(w'*w);
skew_w = [0,-w(3),w(2);...
            w(3),0,-w(1);...
            -w(2),w(1),0];
expSO3 = eye(3) + sin(th)/th*skew_w + (1-cos(th))/th^2*skew_w*skew_w; % Ex. A.11
A = eye(3) + (1-cos(th))/th^2*skew_w + (th-sin(th))/th^3*skew_w*skew_w;
expSE3 = [expSO3,A*v; zeros(1,3),1];    % Ex. A.12

% exp(wt)
rotate_step = expSE3;
wt = [zeros(5,1);phi];
for k = 1:numel(wt)
    rotate_step = subs(rotate_step,se3(k),wt(k));
end

% exp(vt)
prop_step = expSE3;
vt = [0;0;D;Dk;0;0];
for k = 1:numel(wt)
    prop_step = subs(prop_step,se3(k),vt(k));
end

% kinematics
step = prop_step * rotate_step






%% SO(3) Riemannian metric
q1 = sym('q1', [4,1], 'real');
q2 = sym('q2', [4,1], 'real');

% corresponds to  q2 = quat_multiply(dq,q1)
% quat_multiply(q2, quat_inv(q1))
% Space Vehicle Dynamics and Control - Wie
prod_matrix = [q1(4), q1(3), -q1(2), -q1(1); -q1(3), q1(4), q1(1) -q1(2); q1(2), -q1(1), q1(4), -q1(3); q1(1), q1(2), q1(3), q1(4)];
dq = prod_matrix * q2;
th =  2*atan2(norm(dq(1:3)), dq(4));

% 2*atan2((abs(q1x*q2y - q1y*q2x + q1z*q2w - q1w*q2z)^2 + abs(q1x*q2z - q1z*q2x - q1y*q2w + q1w*q2y)^2 + abs(q1x*q2w + q1y*q2z - q1z*q2y - q1w*q2x)^2)^(1/2), q1x*q2x + q1y*q2y + q1z*q2z + q1w*q2w)






%% Convex relaxation of arccos(x) for x\in[-1,1]
% Improved Approximation Algorithms for Maximum Cut and Satisfiability Problems Using Semidefinite Programming - Goemans and Williamson (1995)
x = linspace(-1,1,1e4);
th0 = 2.331122;
alfa = 2/pi/sin(th0);
f1 = acos(x)/(2*pi);
f2 = alfa/4 * (1-x);
f3 = (pi/2-x)/(2*pi);
plot(x,f1,'k', x,f2,'r-.', x, f3, 'b-.');
xlabel('x'); ylabel('arccos(x)/(2*\pi)');
legend('cos^{-1}', 'relaxation', 'line fit');






%% MRP2DCM
% A Survey of Attitude Representations - Malcom D. Shuster (1993)
% Eq. 255b
clear all; clc;
p = sym('p', [3,1], 'real');
skew = [0, -p(3), p(2);...
    p(3), 0, -p(1);...
    -p(2), p(1), 0];
dcm = simplify(eye(3,3) + (8*skew^2 - 4*(1 - norm(p)^2)*skew)/(1+norm(p)^2)^2);






%% Dubins plane
clear all; clc;
syms x y z v psi gamma phi alpha ua upsid ualphad
syms g rho Area mass Cd0 Kd alpha_0
assume(x, 'real');
assume(y, 'real');
assume(z, 'real');
assume(v, 'real');
assume(psi, 'real');
assume(gamma, 'real');
assume(phi, 'real');
assume(alpha, 'real');
assume(ua, 'real');
assume(upsid, 'real');
assume(ualphad, 'real');
assume(g, 'real');
assume(rho, 'real');
assume(Area, 'real');
assume(mass, 'real');
assume(Cd0, 'real');
assume(Kd, 'real');
assume(alpha_0, 'real');

Fl = pi*rho*Area*v^2*alpha;
Fd = rho*Area*v^2*(Cd0 + 4*pi^2*Kd*alpha^2);

X = [x;y;z;psi;v;gamma;phi;alpha];
U = [ua;upsid;ualphad];
Xdot = [v*cos(psi)*cos(gamma);...
        v*sin(psi)*cos(gamma);...
        v*sin(gamma);...
        -Fl*sin(phi)/(mass*v*cos(gamma));...
        ua - Fd/mass - g*sin(gamma);...
        Fl*cos(phi)/(mass*v) - g*cos(gamma)/v;...
        upsid;...
        ualphad];

Adub = jacobian(Xdot, X);
Bdub = jacobian(Xdot, U);
