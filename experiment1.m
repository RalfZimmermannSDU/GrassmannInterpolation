%% Experiment 1: Q factor experiment
function experiment1(n,p,t0,t1,m,LoR,maxsteps)
% n,p      = dimensions of the Grassmann manifold
% t0,t1,m  = time interval and number of points 
% LoR      = apply maxvol on the left- or rightmost endpoint
% maxsteps = maximal number of maxvol-algorithm steps.

rng(123123,'twister')
close all

M = matrix_tools();

% Genrate the curve Y(t)
Y0 = rand(n,p); Y1 = 0.5*rand(n,p); Y2 = 0.5*rand(n,p); Y3 = 0.25*rand(n,p);

Data = cell(1,2);
Data_P = cell(1,2);

dData = cell(1,2);
dDatahor = cell(1,2);
dData_P = cell(1,2);


[Data{1}, dData{1}] = curve2(t0,Y0,Y1,Y2,Y3);
c1 = norm(eye(p)/Data{1}(1:p,1:p) , 'fro');
[Data{2},dData{2}] = curve2(t1,Y0,Y1,Y2,Y3);
c2 = norm(eye(p)/Data{2}(1:p,1:p) , 'fro');

% Obtain horizontal lift to St(n,p)
dDatahor{1} = (dData{1}*Data{1}'+Data{1}*dData{1}')*Data{1};
dDatahor{2} = (dData{2}*Data{2}'+Data{2}*dData{2}')*Data{2};


% Apply maxvol method
if LoR == "L"
    [~,P] = maxvol(Data{1},maxsteps);
    Data_P{1} = P*Data{1};
    Data_P{2} = P*Data{2};
    dData_P{1} = P*dData{1};
    dData_P{2} = P*dData{2};
else
    [~,P] = maxvol(Data{2},maxsteps);
    Data_P{1} = P*Data{1};
    Data_P{2} = P*Data{2};
    dData_P{1} = P*dData{1};
    dData_P{2} = P*dData{2};
end

p_c1 = norm(eye(p)/Data_P{1}(1:p,1:p) , 'fro');
p_c2 = norm(eye(p)/Data_P{2}(1:p,1:p) , 'fro');


% Lagrange interpolation

e1s = [];
e2s = [];
e3s = [];

fe1s = [];
fe2s = [];
fe3s = [];
for i = 1:m+1
    t = (i-1)/m*t1+t0;
    U = curve2(t,Y0,Y1,Y2,Y3);
    Q1 = P'*Interpolate_Gr([t0 t1],Data_P,t,'local_lag'); % MV coords
    Q2 = Interpolate_Gr([t0 t1],Data,t,'local_lag'); % Standard local coords
    Q3 = Interpolate_Gr([t0 t1],Data,t,'normal_lag'); % Normal
    
    NU = norm(U*U','fro');
    e1 = norm(Q1*Q1'-U*U','fro')/NU; % MV coords
    e2 = norm(Q2*Q2'-U*U','fro')/NU; % Standard local coords
    e3 = norm(Q3*Q3'-U*U','fro')/NU; % Normal

    e1s(i) = e1;
    e2s(i) = e2;
    e3s(i) = e3;


    fe1s(i) = norm(Q1'*Q1 - eye(p),'fro'); % MV coords
    fe2s(i) = norm(Q2'*Q2 - eye(p),'fro'); % Standard local coords
    fe3s(i) = norm(Q3'*Q3 - eye(p),'fro'); % Normal

end


ts = linspace(t0,t1,m+1);

f = figure;
f.Position = [40,800,1200*5/6,650*5/6*2];

subplot(2,2,1)
plot(ts,e1s,'-')
hold on
plot(ts,e2s,'--')
plot(ts,e3s,'-.')
legend("MV coords","Local coords","Normal coords")
title("Error (Lagrange)")

xlabel("t")
ylabel("Rel. error")

subplot(2,2,3)
semilogy(ts,fe1s,'-')
hold on
semilogy(ts,fe2s,'--')
semilogy(ts,fe3s,'-.')
xlabel("t")
ylabel("Feasibility")
legend("MV coords","Local coords","Normal coords")
title("Feasibility (Lagrange)")


% Hermite interpolation

e1s = [];
e2s = [];
e3s = [];

fe1s = [];
fe2s = [];
fe3s = [];

for i = 1:(m+1)
    t = (i-1)/m*t1+t0;
    U = curve2(t,Y0,Y1,Y2,Y3);
    Q1 = P'*Interpolate_Gr([t0 t1],Data_P,t,'local_herm',dData_P); % MV coords
    Q2 = Interpolate_Gr([t0 t1],Data,t,'local_herm',dData); % standard local coords
    Q3 = Interpolate_Gr([t0 t1],Data,t,'normal_herm',dDatahor); % Normal 
    
    NU = norm(U*U','fro');
    e1 = norm(Q1*Q1'-U*U','fro')/NU; % MV coords
    e2 = norm(Q2*Q2'-U*U','fro')/NU; % standard local coords
    e3 = norm(Q3*Q3'-U*U','fro')/NU; % Normal 

    e1s(i) = e1;
    e2s(i) = e2;
    e3s(i) = e3;


    fe1s(i) = norm(Q1'*Q1 - eye(p),'fro'); % MV coords
    fe2s(i) = norm(Q2'*Q2 - eye(p),'fro'); % standard local coords
    fe3s(i) = norm(Q3'*Q3 - eye(p),'fro'); % Normal 

end


subplot(2,2,2)
plot(ts,e1s,'-')
hold on
plot(ts,e2s,'--')
hold on
plot(ts,e3s,'-.')
xlabel("t")
ylabel("Rel. error")

legend("MV coords","Local coords","Normal coords")
title("Error (Hermite)")

subplot(2,2,4)
semilogy(ts,fe1s,'-')
hold on
semilogy(ts,fe2s,'--')
semilogy(ts,fe3s,'-.')
legend("MV coords","Local coords","Normal coords")
title("Feasibility (Hermite)")
xlabel("t")
ylabel("Feasibility")
sgtitle("Relative interpolation errors and feasibilities")

fontsize(f,18,"pixels")
exportgraphics(f,"experiment_1.png","Resolution",300);


% Table of norms of inverted upper p x p blocks before and after maxvol
Point = [0,1]';
Before = [c1,c2]';
After = [p_c1,p_c2]';

T = table(Point,Before,After);
disp(T);

end
function [Q,dQ] = curve2(t,Y0,Y1,Y2,Y3)
    M = matrix_tools();
    Y = Y0 + t*Y1 + t^2*Y2 + t^3*Y3;
    dY = Y1 + 2*t*Y2 + 3*t^2*Y3;
    [Q,~] = qr(Y,'econ');
    if nargout < 2
        return
    else
        [dQ,~] = M.dQR(Y,dY);
        return
    end
end