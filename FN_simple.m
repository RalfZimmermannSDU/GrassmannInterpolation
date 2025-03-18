% FH system without diffusion
clear all
close all
% Play with the parameters in the 2D FN system.
% This could be the basis of an exercise in MONA or ODE

% Parametes (not to be compared 1-1 to the FN system with diffusion)

% type II
a = 0.1;
b = 1.5;
epsilon = 0.01;

% Relaxation-Oscillaions
a = 0;
b = 1.45;
epsilon = 0.01;

% System:
% ut = -u^3+u-v
% vt = epsilon(u-bv+a)

% write yt = [ut,vt]'

% Initial value? Try y0 = [0,0]'
 
y = [-0.15,-0.2]';

nt = 10e5;
T = 900;
dt = T / nt;

tvals = 0:dt:T;

Y = [y];
ts = [0];
% Iterate via simple Euler
for i = 1:nt
    t = i*dt;
    
    y(1) = y(1) + dt * (-y(1)^3+y(1)-y(2));
    y(2) = y(2) + dt * (epsilon*(y(1)-b*y(2)+a));

    if mod(i,500) == 0
        Y(:,end + 1) = y;
        ts(end + 1) = t;
    end
end

figure
plot(ts,Y(1,:))

figure
plot(Y(1,:),Y(2,:))