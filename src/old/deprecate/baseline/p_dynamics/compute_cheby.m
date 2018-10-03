function [T, T_dot] = compute_cheby(K,N,t)

T = zeros(N+1,K+1); % up to order N
U = zeros(N+1,K+1); 
T_dot = zeros(N+1,K+1);

T(1,:) = 1.0; %N = 0
T(2,:) = t; %N = 1

U(1,:) = 1.0; %N = 0
U(2,:) = 2*t; %N = 1

T_dot(2,:) = 1.0; %N = 1

% V = c(1)*T(:,1)+c(2)*T(:,2);
% V_dot = c(1)*T_dot(:,1) + c(2)*T_dot(:,2);

for n = 2:N %order, posn = n+1
    T(n+1,:) = 2.0*t.*T(n,:) - T(n-1,:);
    U(n+1,:) = 2.0*t.*U(n,:) - U(n-1,:);
    T_dot(n+1,:) = n*U(n,:);    
end

% V = C*T;
% V_dot = C*T_dot;


end