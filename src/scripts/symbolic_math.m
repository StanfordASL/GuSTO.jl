%% Discretization of double integrator dynamics for ZOH
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






%% Rigid body dynamics on SE(3)
%% A Survey of Attitude Representations - Malcom D. Shuster (1993)
clear all; clc;
syms Jxx Jyy Jzz Jxy Jxz Jyz
assume(Jxx, 'real');
assume(Jyy, 'real');
assume(Jzz, 'real');
p = sym('p', [3,1], 'real');
wb = sym('wb', [3,1], 'real');
M = sym('M', [3,1], 'real');
J = diag([Jxx,Jyy,Jzz]);

pb_dot = 0.25*((1-p'*p)*wb + 2*cross(p,wb) + 2*dot(wb,p)*p);  % Eq. 338
wbdot = inv(J)*(M - cross(wb,J*wb));
Apb = [jacobian(pb_dot, [p;wb]); ...
    jacobian(wbdot, [p;wb])];
Bpb = [jacobian(pb_dot, [M]); ...
    jacobian(wbdot, [M])];






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
