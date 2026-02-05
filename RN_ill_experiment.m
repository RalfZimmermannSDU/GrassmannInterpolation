%% Experiment: Ill behaved RN coordinates

clear all
close all
clc;

rng default

n = 10;
p = 3;
%x = linspace(0,1,180000);

M = matrix_tools();

U0 = [eye(p); zeros(n-p,p)];
Dir = M.vectorG(U0);
v = M.vectorG(U0)*0.5;

Y0 = rand(n,p); Y1 = rand(n,p); Y2 = rand(n,p);
c = @(t) U0 + Y0*t + t^2*Y1/2 + t^3*Y2/3;
dc = @(t) Y0 + t*Y1 + t^2*Y2;

q = 1/sqrt(2) * [eye(p);eye(p);zeros(n-2*p,p)];


t = 0.5;


%[n,p] = size(Y);

%[Q,R] = qr(Y,'econ');

ts = linspace(1e-13,eps,100);
%ts = linspace(1e-11,eps,100);
norms_RN = [];
norms_loc = [];
for i = 1:100
    t = ts(i);
    %eps = 0;

    [Q,R] = qr(c(t),'econ');
    
    D = diag(sign(diag(R)));
    D(D == 0) = 1;
    
    Q = -Q;
    R = -R;

    dY = dc(t);

    PL = tril(ones(p,p),-1);
    L = PL.*(Q'*dY / R);
    X = L - L';
    
    dR = Q'*dY - X*R;
    dQ = (eye(n) - Q*Q') * dY / R + Q*X;


    %U = c(eps);
    %dU = dc(eps);
    U = Q;
    dU = dQ;
    dUc = (dU*U' + U*dU')*U; % Hor. lift

    % Derivative of Log:
    B = (eye(n) - q*q') * U /(q'*U);
    dB = (eye(n) - q*q') *(dUc - U * ((q'*U) \ (q'*dUc) ))/(q'*U);

    [U,S,V] = svd(B,'econ');
    [dU,dS,dV] = M.dSVD(U,S,V,dB);

    dLog = dU * atan(S) * V' + U * (1./(1+S.^2).* dS)  * V' + U * atan(S) * dV';

    % FD approx
    h = 0.000000001;
    lmh = GrassLog(q,GrassExp(Q,-h*dUc,1));
    lph = GrassLog(q,GrassExp(Q,h*dUc,1));
    FD = (lph - lmh)/(2*h);
    normDiff = norm(FD - dLog,'fro') / norm(FD,'fro');
    norms_RN(i) = normDiff;

    % The local coordiante version
    % q is the center
    [UQ,~,~] = M.house_qr(q);
    [W,Y] = M.house_block(UQ);
    q_centered = M.applyWYT(q,W,Y);
    Q_centered = M.applyWYT(Q,W,Y);
    dQc_centered = M.applyWYT(dUc,W,Y);

    dQ_local = M.dLocalCoordG(Q_centered,dQc_centered,n,p);

    Bl = M.LocalCoordG(M.applyWYT(GrassExp(Q,-h*dUc,1),W,Y),n,p);
    Bh = M.LocalCoordG(M.applyWYT(GrassExp(Q,h*dUc,1),W,Y),n,p);
    normDiff = norm((Bh-Bl)/(2*h)-dQ_local,'fro') / norm((Bh-Bl)/(2*h),'fro');
    norms_loc(i) = normDiff;
end
%%
f = figure;
f.Position = [40,800,1200*6/5*1/2,650*5/5];
set(f, 'DefaultTextInterpreter', 'latex')
lw = 3;
semilogy(ts,norms_loc,':','LineWidth',lw,'Color','#1171BE')
hold on
semilogy(ts,norms_RN,'LineWidth',lw,'Color','#DD5400')

for i = 1:100
    if isnan(norms_RN(i))
        semilogy(ts(i), 1e-1, 'o','LineWidth',4,'Color','#DD5400'); % Add a star for NA values
    end
end

ylim([1e-9,10])

title("Relative error of FD approximations")
legend("CCCs","RN")
xlabel("t")
ylabel("Rel. error")

fontsize(20,"points")

exportgraphics(f,"fig_rn_log.png","Resolution",600)
% [U,R] = qr(c(eps),'econ');
% 
% 
% h = 0.00000001;
% eps = eps + h;
% [Uh,Rh] = qr(c(eps),'econ');


%norm((Uh - U)/h - dQ)
function L = GrassLog(q,p)
% This is the classic logarithm 
[n,~] = size(p);
mat = (eye(n) - q*q') * (p /(q'*p)); % inefficient, but readable

[U,S,V] = svd(mat,'econ');

L = U * atan(S) * V';
end
function Y = GrassExp(U,dir,t)
    % This is equation 3.10 in Bendokat and Zimmermann 2024
    if nargin < 3 || isempty(t)
        t = 1;
    end
    Dir = t * dir;
    [Q,S,V] = svd(Dir,'econ');
    Scos = diag(cos(diag(S)));
    Ssin = diag(sin(diag(S)));

    Y = U * (V*Scos * V') + Q * (Ssin * V');

    %[Y,~] = qr(Y,'econ');
end

function [C] = vec2skew(c_vec)
% build skew-symmetric matrix C from c_vec
% recall that q = dim(c_vec) = 0.5(p*(p-1))
    q = length(c_vec);
    p = floor(0.5 + sqrt(0.25 + 2*q));
    C = zeros(p,p);
    count = 1;
    for j = 1:p
        for k = j+1:p
            C(j,k) = c_vec(count);
            C(k,j) = -C(j,k);
            count = count+1;
        end
    end
    return;
end