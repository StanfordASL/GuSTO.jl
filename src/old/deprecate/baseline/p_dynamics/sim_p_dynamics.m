clear all; close all; clc;
warning('off','MATLAB:nargchk:deprecated');
addpath('../quat_functions');
         
load_sim_config;

% uses Tomlab to set up LQR style PSM set-up 
[MP_Prob,L_e_mp,MP_st] = setup_MP(n,m,...
    f,B,df,...
    state_constr,ctrl_constr,c_nl,...
    N_mp,Tf,dt,...
    x_eq,u_eq,P,Q,R,...
    xnorm,dxnorm,'MP');

if exist(filename, 'file') == 2
    load(filename);
else
    mp_warm = struct('Tp',Tf,'shift',0,'sol',0,...
                's_t',MP_st,'state',[],'ctrl',[],'result',[]);
end
mp_warm = struct('Tp',Tf,'shift',0,'sol',0,...
            's_t',MP_st,'state',[],'ctrl',[],'result',[]);

tic
[MP_state,MP_ctrl,converged_MP,mp_warm] = compute_MP(MP_Prob,...
    test_state,test_state,state_constr,ctrl_constr,x_eq,u_eq,...
    n,m,N_mp,L_e_mp,mp_warm);
toc
disp('MP:'); disp(converged_MP);

mp_warm.sol = 1;
save(filename,'mp_warm');
