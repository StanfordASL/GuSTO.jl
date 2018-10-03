clear all
close all
clc

% Converts from Euler-Rodrigues Symmetric Parameters to angular velocities
% and so on

N = 100;
T = 5; % horizon, in [s]
dt = T/N;
DIM_ROT = 3;

%% Generate samples
p = zeros(DIM_ROT, N);

A = zeros(N,N);
A = A - 30.0/12.0*diag(ones(N,1));
A = A + 16.0/12.0*diag(ones(N-1,1),1) + 16.0/12.0*diag(ones(N-1,1),-1);
A = A - 1/12.0*diag(ones(N-2,1),2) - 1/12.0*diag(ones(N-2,1),-2);
A = A./(dt^2);


Rinv = inv(A'*A);
% Rinv = (Rinv+ Rinv')/2; % Make it symmetrical, it's not before (rounding errors)

M = Rinv;
for col_i = 1:size(Rinv,2)
    max_col_i = max(Rinv(:, col_i));
    M(:,col_i) = Rinv(:,col_i) / (max_col_i*N);
end

figure(1)
subplot(2,1,1);
for col_i = 1:size(Rinv,2)
    grid on; hold on;
    plot(Rinv(:,col_i));
end
subplot(2,1,2);
for dim=1:DIM_ROT
    eps = randn(1,N)*sqrtm(Rinv);
    Qperturbed = eps;
    p(dim,:) = eps;
    grid on; hold on;
    plot(Qperturbed);
end
movegui('west')



%% Compute Euler rates and accelerations using finite differencing
p_dot  = zeros(DIM_ROT,N);
p_ddot = zeros(DIM_ROT,N);
for dim=1:DIM_ROT
    p_dot(dim,:) = ([diag(ones(N,1))-diag(ones(N-1,1),-1)] * p(dim,:)')';
    p_ddot(dim,:) = (A * p(dim,:)')';
end

figure(2)
for derivative_id = 1:3
    subplot(3,1,derivative_id); grid on; hold on;
    if derivative_id ==1
        vec = p;
    elseif derivative_id ==2
        vec = p_dot;
    else
        vec = p_ddot;
    end
    for dim=1:DIM_ROT
        plot(vec(dim,:));
    end
end


%% Compute angular velocity, angular acceleration
clc
omega = zeros(DIM_ROT, N);
omega_dot = zeros(DIM_ROT, N);
for nidx = 1:N
    pn = p(:,nidx);
    omega(:,nidx) = (4/(1+sum(pn.^2))) * (eye(3)+brackets(pn))/(eye(3)-brackets(pn)) * p_dot(:,nidx);
end
omega_dot = [zeros(3,1), omega(2:N)-omega(:,1:N-1)];
omega_dot;

%% Compute torques
J_R = ones(3,3);
M = J_R * omega_dot + skew(omega) * J_R * omega;
M

function [mat] = brackets(u)
    mat = [0, u(3), -u(2); ...
           -u(3), 0, u(1); ...
           u(2), -u(1), 0];
end
function S = skew(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end