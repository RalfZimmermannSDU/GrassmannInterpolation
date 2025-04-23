% Toy examples
%
% April 2025
%


clear
rng(123123,'twister')
close all
M = matrix_tools();
% Genreal experiment parameters
lt0 = 0; lt1 = 10; % Time interval for Lagrange wave plot.
t0 = 0; t1 = 1; % Time interval for Hermite interpolation.

n = 1000;
p = 10;

Y0 = rand(n,p); Y1 = 0.5*rand(n,p); Y2 = 0.5*rand(n,p); Y3 = 0.25*rand(n,p);
Data = cell(1,2);
Data_P = cell(1,2);
Data{1} = curve2(t0,Y0,Y1,Y2,Y3);
cond(Data{1}(1:p,1:p))
Data{2} = curve2(t1,Y0,Y1,Y2,Y3);
cond(Data{2}(1:p,1:p))

[~,P] = maxvol(Data{1});
Data_P{1} = P*Data{1};
Data_P{2} = P*Data{2};
cond(Data_P{1}(1:p,1:p))
cond(Data_P{2}(1:p,1:p))

t = 0.5;
Q1 = P'*Interpolate_Gr([t0 t1],Data_P,t,'local_lag');
U = curve2(t,Y0,Y1,Y2,Y3);
%c1 = cond(U(1:p,1:p));
e1 = norm(Q1*Q1'-U*U','fro')

Q2 = Interpolate_Gr([t0 t1],Data,t,'local_lag');
U = curve2(t,Y0,Y1,Y2,Y3);
%c2 = cond(U(1:p,1:p))
e2 = norm(Q2*Q2'-U*U','fro')


e1s = [];
e2s = [];
for i = 1:101
    t = (i-1)/100;
    U = curve2(t,Y0,Y1,Y2,Y3);
    Q1 = P'*Interpolate_Gr([t0 t1],Data_P,t,'local_lag'); % Well-behaved
    Q2 = Interpolate_Gr([t0 t1],Data,t,'local_lag'); % Less well-behaved
    e1 = norm(Q1*Q1'-U*U','fro');
    e2 = norm(Q2*Q2'-U*U','fro');

    e1s(i) = e1;
    e2s(i) = e2;

end
figure
plot(0:0.01:1,e1s)
hold on
plot(0:0.01:1,e2s)

legend("Well-behaved","Less well-behaved")

Data = cell(1,2);
dData = cell(1,2);
Data_P = cell(1,2);
dData_P = cell(1,2);

[Data{1}, dData{1}] = curve2(t0,Y0,Y1,Y2,Y3);
cond(Data{1}(1:p,1:p))
[Data{2},dData{2}] = curve2(t1,Y0,Y1,Y2,Y3);
cond(Data{2}(1:p,1:p))

[~,P] = maxvol(Data{1});
Data_P{1} = P*Data{1};
Data_P{2} = P*Data{2};
dData_P{1} = P*dData{1};
dData_P{2} = P*dData{2};
cond(Data_P{1}(1:p,1:p))
cond(Data_P{2}(1:p,1:p))

for i = 1:101
    t = (i-1)/100;
    U = curve2(t,Y0,Y1,Y2,Y3);
    Q1 = P'*Interpolate_Gr([t0 t1],Data_P,t,'local_herm',dData_P); % Well-behaved
    Q2 = Interpolate_Gr([t0 t1],Data,t,'local_herm',dData); % Less well-behaved
    e1 = norm(Q1*Q1'-U*U','fro');
    e2 = norm(Q2*Q2'-U*U','fro');

    e1s(i) = e1;
    e2s(i) = e2;

end
figure
plot(0:0.01:1,e1s)
hold on
plot(0:0.01:1,e2s)

legend("Well-behaved","Less well-behaved")

% norm(M.LocalCoordG(Data{1},n,p)-M.LocalCoordG(Data{2},n,p),'fro')
% M.LC_distbound(Data{1},Data{2})
% 
% norm(M.LocalCoordG(Data_P{1},n,p)-M.LocalCoordG(Data_P{2},n,p),'fro')
% M.LC_distbound(Data_P{1},Data_P{2})

% dist_true = [];
% dist_bound = [];
% for i = 1:100
%     U = M.RandS(n,p);
%     V = M.RandS(n,p);
%     U1 = U(1:p,1:p);
%     V1 = V(1:p,1:p);
% 
%     dist_true(i) = norm(M.LocalCoordG(U,n,p)-M.LocalCoordG(V,n,p),'fro')
%     dist_bound(i) = M.LC_distbound(U,V)
% end
% [dist_true,I] = sort(dist_true,'ascend');
% figure
% semilogy(dist_true)
% hold on
% semilogy(dist_bound(I))
% legend("True distance","Bound")







% % Ill--conditioning?
% % Change the left point slightly to the left by a factor h = 0.01 and
% % observe if the the results dont change too much
% h = 0.01;
% Data2 = cell(1,2);
% Data_P2 = cell(1,2);
% Data2{1} = curve2(0-h,Y0,Y1,Y2,Y3);
% cond(Data2{1}(1:p,1:p))
% Data2{2} = curve2(1,Y0,Y1,Y2,Y3);
% cond(Data2{2}(1:p,1:p))
% 
% t0 = t0 - h;
% t1 = t1; 
% 
% [~,P2] = maxvol(Data2{1});
% Data_P2{1} = P*Data2{1};
% Data_P2{2} = P*Data2{2};
% cond(Data_P2{1}(1:p,1:p))
% cond(Data_P2{2}(1:p,1:p))
% 
% % Data and Data2 distance
% d1 = norm(Data{1}*Data{1}'-Data{2}*Data{2}')
% d2 = norm(Data2{1}*Data2{1}'-Data2{2}*Data2{2}')
% 
% t = 0.5;
% 
% U = curve2(t,Y0,Y1,Y2,Y3);
% Q11 = P'*Interpolate_Gr([t0+h t1],Data_P,t,'local_lag'); % Well-behaved
% Q21 = Interpolate_Gr([t0+h t1],Data,t,'local_lag'); % Less well-behaved
% e11 = norm(Q11*Q11'-U*U','fro')
% e21 = norm(Q21*Q21'-U*U','fro')
% 
% Q12 = P2'*Interpolate_Gr([t0 t1],Data_P2,t,'local_lag'); % Well-behaved
% Q22 = Interpolate_Gr([t0 t1],Data2,t,'local_lag'); % Less well-behaved
% e12 = norm(Q12*Q12'-U*U','fro')
% e22 = norm(Q22*Q22'-U*U','fro')
% 
% norm(Q12*Q12' - Q11*Q11')
% norm(Q22*Q22' - Q21*Q21')

function U = curve1(t,W)
    p = 15;
    [U,~,~] = svd(W*(t+1));
    %U = U(:,1:p);
    ps = 15;
    [U,~,~] = svd(U(:,1:ps)*U(:,1:ps)');
    %[U,~,~] = svd(U*U');
    U = U(:,1:p);
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
function [U,P] = maxvol(U)
    [n,p] = size(U);
    
    Usquare = U(1:p,1:p);
    cond_start = cond(Usquare);

    E = sparse(eye(n));
    
    E2 = E;
    warning('off','MATLAB:nearlySingularMatrix')
    for k = 1:30
        B = U / Usquare;
        [b,I] = max(abs(B),[],'all');
        if B(I)<0
            b = -b;
        end
        %disp(num2str(b))
        if abs(b) > 1
            [i,j] = find(~(B-ones(n,p)*b));
            %disp("(i,j) = " + num2str(i) + ", " +num2str(j))
            U = U + (E(:,j) - E(:,i))*(U(i,:)-U(j,:));
            
            Ei = E2(i,:);
            Ej = E2(j,:);
    
            E2(i,:) = Ej;
            E2(j,:) = Ei;
        end
        Usquare = U(1:p,1:p);
        if abs(b) < 1 + 10e-3
            break
        end
    end
    cond_end = cond(Usquare);
    P = E2;
    disp("Maxvol algorithm:")
    disp("num. iter " + num2str(k));
    disp("Condition number before " + num2str(cond_start))
    disp("Condition number after  " + num2str(cond_end))

end

% ps = 300;
% n = 500; 
% p = 15; % Manifold dimension.
% 
% D = [];
% for i = 1:ps/3
%     D(3*i) = (rand()*0.15)^2;
% end
% D = diag(D);
% 
% W = rand(n,ps) * D;
% 
% Data = cell(1,2);
% Data{1} = curve1(0,W);
% cond(Data{1} (1:p,1:p))
% Data{2} = curve1(0.1,W);
% cond(Data{2}(1:p,1:p))
% 
% [U,P] = maxvol(Data{1});
% 
% Data{1} = P*Data{1};
% Data{2} = P*Data{2};
% 
% cond(Data{1}(1:p,1:p))
% cond(Data{2}(1:p,1:p))
% t0 = 0;
% t1 = 0.1;
% norm(Data{1}*Data{1}'-Data{2}*Data{2}')
% err = [];
% for i = 1:11
%     t = (i-1)/10*t1;
%     Q = Interpolate_Gr([t0 t1],Data,t,'local_lag');
%     U = curve1(t,W);
%     cond(U(1:p,1:p));
%     err(i) = norm(Q*Q'-U*U','fro');
% end