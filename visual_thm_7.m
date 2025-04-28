clear
close all
rng(123123,'twister')
n = 1000;
p = 20;

M = matrix_tools();

U = M.RandG(n,p);
Delta = M.vectorG(U);
Delta = Delta/(norm(Delta,'fro')*0.5);
[U,P] = maxvol(U);

B = M.LocalCoordG(U,n,p);
cond_d_bound = M.dphi_cond_bound(B);
Delta = P*Delta/(norm(Delta,'fro')*0.5);

man_dist = [];
bound = [];
I = [];
for i = 1:101
    t = i/4001;
    man_dist(i) = t;
    
    V = M.ExpG(U,Delta,t);
    Btilde = M.LocalCoordG(V,n,p);
    bound(i) = asin(cond_d_bound*norm(B-Btilde,'fro'));

    if ~(t < bound(i))
        disp("Error");
    end

end

% figure
% plot(man_dist,bound)

m = 2000;

U = M.RandG(n,p);
[U, P] = maxvol(U);
B = M.LocalCoordG(U,n,p);
cond_d_bound = M.dphi_cond_bound(B);
Data = zeros(m,3);
conds = [];
for i = 1:m
    Delta = M.vectorG(U);
    Delta = P*Delta / (norm(Delta,'fro')*0.5);
    t = rand()*0.025;

    V = M.ExpG(U,Delta,t);
    Btilde = M.LocalCoordG(V,n,p);
    conds(i) = cond(V(1:p,1:p),'fro');
    Data(i,1) = t;
    Data(i,2) = asin(cond_d_bound*norm(B-Btilde,'fro')/2);
    Data(i,3) = norm(B-Btilde,'fro')/2;
    if ~(t < Data(i,3))
        disp("Error");
    end
end

f = figure;
f.Position = [40,800,1200*5/6,650*5/6];

subplot(1,2,1)
plot(Data(:,1),Data(:,3),'x')
grid on
hold on
xlabel("Manifold distance")
ylabel("Distance between local coords")

title("dist(U,V) vs.  ||B_1-B_2||_0")

subplot(1,2,2)
plot(Data(:,1),Data(:,2),'x')
grid on
xlabel("Manifold distance")
ylabel("Bound")
title("dist(U,V) vs. arcsin(|C ||B_1-B_2||_0)")
fontsize(f,15,"pixels")

exportgraphics(f,"theorem_7.png","Resolution",300);

% figure
% plot(Data(:,1),conds,'x')

function [U, P] = maxvol(U)
    [n,p] = size(U);
    
    Usquare = U(1:p,1:p);
    cond_start = cond(Usquare,'fro');

    E = sparse(eye(n));
    
    E2 = E;
    warning('off','MATLAB:nearlySingularMatrix')
    for k = 1:40
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
    cond_end = cond(Usquare,'fro');
    P = E2;
    disp("Maxvol algorithm:")
    disp("num. iter " + num2str(k));
    disp("Condition number before " + num2str(cond_start))
    disp("Condition number after  " + num2str(cond_end))

end