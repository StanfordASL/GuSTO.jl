function [MP_Prob,L_e_full,s_t] = ...
    setup_MP(n,m,m_f,...
             f,B,df,...
             state_con,u_con,f_con,...
             A,l_con,nl_con,...
             N,Tf,dt,...
             x_eq,u_eq,f_eq,...
             P,Q,R,G,...
             ts_wps, wps,...
             qnorm,dqnorm,Name)

%% Constants

%State bounds
x_L = state_con(:,1);
x_U = state_con(:,2);

%Control bounds
u_L = u_con(:,1);
u_U = u_con(:,2);

f_L = f_con(:,1);
f_U = f_con(:,2);

%Number of collocation points
K = N;

%CGL nodes
[s_t,w] = clencurt(K); %t_t: [-1, 1] : <-> : [0, Tf]
s = fliplr(s_t); %t: [1, -1]

%% Final solution interpolation matrix

tau_full = 0:dt:Tf; 
s_e_full = (2*tau_full - Tf)/Tf; %[-1, 1]

idxs = zeros(numel(ts_wps),1);
for k = 1:numel(ts_wps)
  t_wp = ts_wps(k);
  [val,idxs(k)] = min(abs(s_t-t_wp)); % TODO(acauligi): check if s_t or s here
end

%Lagrange polynomial evaluation at the interpolation points
% L_e = compute_Lagrange(length(s_e)-1,N,s_e,s_t);
L_e_full = compute_Lagrange(length(s_e_full)-1,N,s_e_full,s_t);

%% Get Differentiation matrix

D = ChebyshevDiffMatrix(N,s); %arranged for forward time
D = kron(D,eye(n));

%% Variables

%State node values: [x_t0,...,x_tN]
%Control node values: [u_t0,...,u_tN]

% n_vars = (N+1)*(n+m);

%% Define problem

% u_eq = zeros(m,1);

x_eq_all = kron(ones(N+1,1),x_eq);
u_eq_all = kron(ones(N+1,1),u_eq);
f_eq_all = kron(ones(N+1,1),f_eq);

xu_eq = [x_eq_all;u_eq_all;f_eq];

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
    
for k = 1:numel(idxs)
  idx = idxs(k);
  xu_L(n*(idx-1)+1:n*(idx-1)+4) = wps(1:4,k);
  xu_U(n*(idx-1)+1:n*(idx-1)+4) = wps(1:4,k);
end

% linear equality constraints
b_L = l_con(:,1);
b_U = l_con(:,2);

% nonlinear equality constraints
c_L = nl_con(:,1);
c_U = nl_con(:,2); 
       
MPC_cost = @(xu) (Tf/2)*(xu-xu_eq)'*F*(xu-xu_eq);
MPC_grad = @(xu) Tf*F*(xu-xu_eq);
MPC_hess = @(xu) Tf*F;

xu0 = zeros((n+m+m_f)*(N+1),1);

% Name = 'MP';
MP_Prob = conAssign(MPC_cost,MPC_grad,MPC_hess,[],...
            xu_L,xu_U,Name, xu0,...
            [], 0, A, b_L, b_U,...
            'MP_con','MP_conJ',[],[],...
            c_L,c_U,...
            [],[],[],[]);
        
MP_Prob.user.x_act = zeros(n,1);
MP_Prob.user.D = D;
MP_Prob.user.n = n;
MP_Prob.user.m = m;
MP_Prob.user.N = N;
MP_Prob.user.f = f;
MP_Prob.user.df = df;
MP_Prob.user.B_full = B_full;
MP_Prob.user.qnorm = qnorm;
MP_Prob.user.dqnorm = dqnorm;
MP_Prob.user.Tf = Tf;
MP_Prob.user.P = P;
MP_Prob.user.x_eq = x_eq;

end
