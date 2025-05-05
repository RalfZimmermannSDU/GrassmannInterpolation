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

% Flags
snapshots_FN = 1;

% Experiments to run
run_exp_1 = 1;
run_exp_2 = 1;


%% Experiment 1: Q factor experiment
n = 1000;
p = 10;
t0 = 0; t1 = 1;
m = 100;
LoR = 'R';
maxsteps = 30;
if run_exp_1
    experiment1(n,p,t0,t1,m,'R',maxsteps);
end