clear all; clc; close all;

% Just some notes, since the spy command is easier to use here
N = 100;
T = 1; % horizon, in [s]
dt = T/N;
mid_N = floor(N/2); N_1sthalf = mid_N; N_2ndhalf = N-mid_N;

%% combinatioon of forwards/backwards 2nd order difference'
disp('How it should be:')
A = +diag(ones(N,1))-2*diag(ones(N-1,1),1)+diag(ones(N-2,1),2);
A(end-1:end, end-2:end) = [[1;1],[-2;-2],[1;1]];
A
det(A)
det(A'*A)
%% or (as before, but with middle point)
A(1:mid_N,1:mid_N) = diag(ones(N_1sthalf,1)) -2*diag(ones(N_1sthalf-1,1),1) +diag(ones(N_1sthalf-2,1),2);
A(mid_N+1:end,mid_N+1:end) = diag(ones(N_2ndhalf,1)) -2*diag(ones(N_2ndhalf-1,1),-1) +diag(ones(N_2ndhalf-2,1),-2);
A(mid_N,:)=zeros(size(A(mid_N,:)));
A(mid_N,mid_N-1:mid_N+1) = [1,-2,1]; % make a junction in the middle
A(mid_N+1,mid_N:mid_N+2) = [1,-2,1];
A(mid_N+2,mid_N:mid_N+2) = [1,-2,1];
A

%%  or (from paper but is wrong!!!)
A(1:mid_N,1:mid_N) = diag(ones(N_1sthalf,1)) -2*diag(ones(N_1sthalf-1,1),-1) +diag(ones(N_1sthalf-2,1),-2);
A(mid_N+1:end,mid_N+1:end) = diag(ones(N_2ndhalf,1)) -2*diag(ones(N_2ndhalf-1,1),1) +diag(ones(N_2ndhalf-2,1),2);
% A(mid_N,:)=zeros(size(A(mid_N,:)));
% A(mid_N,mid_N-1:mid_N+1) = [1,-2,1]; % make a junction in the middle
% A(5,6) = 0; A(6,5)=0;
% A(4,6) = 0; A(5,7)=0;

% A = 2*diag(ones(N,1)) -diag(ones(N-1,1),1)     -diag(ones(N-1,1),-1);
% A =   diag(ones(N,1)) -2*diag(ones(N-1,1),-1)  +diag(ones(N-2,1),-2);
% A = 2*diag(ones(N,1)) - (diag(ones(N-1,1),1)+diag(ones(N-1,1),-1));

R = A'*A;
Rinv = inv(R);
% Rinv = (Rinv+ Rinv')/2; % Make it symmetrical, it's not before (rounding errors)


close all
    subplot(1,2,1);
for col_i = 1:size(Rinv,2)
    grid on; hold on;
    plot(Rinv(:,col_i));
end
    subplot(1,2,2);
for i=1:10
    eps = randn(1,N)*Rinv;
    Qperturbed= eps;
    grid on; hold on;
    plot(Qperturbed);
end
movegui('west')

%% other (CURRENT IMPLEMENTATION!)
A = zeros(N,N);
A = A - 30.0/12.0*diag(ones(N,1));
A = A + 16.0/12.0*diag(ones(N-1,1),1) + 16.0/12.0*diag(ones(N-1,1),-1);
A = A - 1/12.0*diag(ones(N-2,1),2) - 1/12.0*diag(ones(N-2,1),-2);
A = 100*A./(dt^2);



Rinv = inv(A'*A);
% Rinv = (Rinv+ Rinv')/2; % Make it symmetrical, it's not before (rounding errors)

M = Rinv;
for col_i = 1:size(Rinv,2)
    max_col_i = max(Rinv(:, col_i));
    M(:,col_i) = Rinv(:,col_i) / (max_col_i*N);
end



close all
subplot(2,1,1);
for col_i = 1:size(Rinv,2)
    grid on; hold on;
    plot(Rinv(:,col_i));
end
subplot(2,1,2);
for i=1:10
%     eps = mvnrnd(zeros(1,N),Rinv,1)
%     eps = randn(1,N)*Rinv;
    eps = randn(1,N)*sqrtm(Rinv);
    Qperturbed = eps;
    grid on; hold on;
    plot(Qperturbed);
end
movegui('west')

%% as in paper, kind of 
%  for 1st half, do backwards, for 2nd, do forwards
A = zeros(N,N);
% A(1:mid_N,1:mid_N) = diag(ones(N,1)) -2*diag(ones(N-1,1),-1) +diag(ones(N-2,1),-2);
A(1:mid_N,1:mid_N) = diag(ones(N_1sthalf,1)) -2*diag(ones(N_1sthalf-1,1),1) +diag(ones(N_1sthalf-2,1),2);
A(mid_N+1:end,mid_N+1:end) = diag(ones(N_2ndhalf,1)) -2*diag(ones(N_2ndhalf-1,1),-1) +diag(ones(N_2ndhalf-2,1),-2);