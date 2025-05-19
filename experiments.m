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
snapshots_FN = 0;

% Experiments to run
run_exp_1 = 0;
run_exp_2 = 1;


%% Experiment 1: Q factor experiment
n = 1000;
p = 10;
t0 = 0; t1 = 1;
m = 200;
LoR = "R";
maxsteps = 3;
if run_exp_1
    experiment1(n,p,t0,t1,m,LoR,maxsteps);
end


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
    FN_create_snapshots(p,points,num_time_pts)
end

p = 8; % has to be <= p_snap.

% Load data
Data = load("snapshots_FN_model/snapshot_N_91.mat"); % Used as interpolation data
Data_ref = load("snapshots_FN_model/snapshot_N_501.mat"); % Used for high_res plots

m2 = 51;
if run_exp_2
    FN_interpolate(p,Data,Data_ref,points,m2,t0,t1,h,h2,maxsteps)
end
















