%% Numerical experiments 
%
% This script allows the user to run the numerical experiments presented in
% the paper. 
% 
% If you have not generated the POD-data for the
% FitzHugh-Naragumo system, then be aware that doing so takes a while.
% To generate the snapshots, set `snapshots_FN' to 1. This will generate a
% folder where the data subsequently will be placed. 
%
clear 
close all

% Flags
snapshots_FN = 1;

% Experiments to run
run_exp_1 = 1;
run_exp_2 = 0;


%% Experiment 1: Q factor experiment
n = 1000;
p = 10;
t0 = 0; t1 = 1;
m = 200;
LoR = "R";
maxsteps = 35;
if run_exp_1
    [ts, e1s, e2s, e3s, e11s, e22s, e33s] = experiment1(n,p,t0,t1,m,LoR,maxsteps);
end

%% Figure
f = figure;
f.Position = [40,800,1200*6/5,650*5/5];
set(f, 'DefaultTextInterpreter', 'latex')
lw = 3;
subplot(1,2,1)
plot(ts,e1s,'-','LineWidth',lw)
hold on
plot(ts,e2s,'--','LineWidth',lw)
plot(ts,e3s,'-.','LineWidth',lw)
legend("MV coords","Local coords","Normal coords",'Interpreter','latex')
title("Error (Lagrange)",'Interpreter','latex')

xlabel("t",'Interpreter','latex')
ylabel("Rel. error",'Interpreter','latex')

subplot(1,2,2)
plot(ts,e11s,'-','LineWidth',lw)
hold on
plot(ts,e22s,'--','LineWidth',lw)
hold on
plot(ts,e33s,'-.','LineWidth',lw)
xlabel("t",'Interpreter','latex')
ylabel("Rel. error",'Interpreter','latex')

legend("MV coords","Local coords","Normal coords",'Interpreter','latex')
title("Error (Hermite)",'Interpreter','latex')
fontsize(18,"points")
%%
exportgraphics(f,"fig_4.png","Resolution",600)

%% Experiment 2: FN system

t0 = 0.03;
t1 = 0.08;
h = 0.01;
h2 = 0.0001;
m = 101;
points = linspace(t0,t1,m);

% n should not be changed, but p can be varied to a small extend. 
n = 1024;
p_snap = 12; % Number of left singular vectors to store
num_time_pts = 10e5; % number of total time steps. 

% The snapshot data is also used to compute interpolation errors, to it
% makes sense to compute at least some for testing.
if snapshots_FN
    points = t0:h:t1;
    FN_create_snapshots(points,num_time_pts)
end
%%
p = 8; % has to be <= p_snap.

% Load data
Data = load("snapshots_FN_model/snapshot_N_91.mat"); % Used as interpolation data
Data_ref = load("snapshots_FN_model/snapshot_N_501.mat"); % Used for high_res plots

m2 = 51;
if run_exp_2
    FN_interpolate(p,Data,Data_ref,points,m2,t0,t1,h,h2,maxsteps)
end

%%
clear all
close all
M = matrix_tools();
Data = load("snapshots_FN_model/snapshot_N_6.mat")
Data_ref = load("snapshots_FN_model/snapshot_N_91.mat")
Data_ref = Data_ref.data_u(1:51);
mh = 10;
Udata = cell(1,6);
dUdata = cell(1,6);

nx = 1024;
% Prepare data
psel = 8;
for i = 1:6
    [U,S,V] = svd(Data.Data{i}(1:nx,:),'econ');
    [k,~] = size(S);
    S = diag(S);
    tol = 10e-8;
    p = sum(S>tol);
    fprintf("Need to truncate to %2d first columns\n",p)
    S = diag(S(1:p));
    [dU,~, ~] = M.dSVD(U(:,1:p),S,V,Data.Data_dot{i}(1:nx,:));
    
    U = U(:,1:psel);
    dU = dU(:,1:psel);

    dU = (dU*U'+U*dU')*U;
    Udata{i} = U;
    dUdata{i} = dU;
end
%%
for i = 1:6
    [U,S,V] = svd(Data.Data{i}(1:nx,:),'econ');
    [Uh,Sh,Vh] = svd(Data.Data_ph{i}(1:nx,:),'econ');

    h = 0.0000001;
    U = U(:,1:psel);
    Uh = Uh(:,1:psel);
    dUFD = (Uh - U) / h;
    dUFD = (dUFD*U'+U*dUFD')*U;
    norm(dUFD - dUdata{i},'fro') 
    %dUdata{i} = dUFD;
end

%%
Data_ref = load("snapshots_FN_model/snapshot_N_501.mat")
Data_ref = Data_ref.data_u(1:501);
mh = 100;
maxsteps = 30;
FN_interpolate_update(Udata, dUdata, Data_ref, 0.03:0.01:0.08, mh, psel, maxsteps);