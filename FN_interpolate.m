% Interpolate the FN system
clear 
close all

Data = load("snapshots_FN_model/snapshot_N_91.mat");


% First, consider two points
points = 0.03:0.001:0.12;
Ias = [0.03 0.05];

p = 8;

I = ismember(points,Ias);

n = sum(I);
[~,m] = size(points);

% Extract data
u_data = cell(1,n);
v_data = cell(1,n);
u_dot_data = cell(1,n);
v_dot_data = cell(1,n);

k = 1;
for i = 1:m
    if I(i)
        u_data{k} = Data.data_u{i}(:,1:p);
        v_data{k} = Data.data_v{i}(:,1:p);
        u_dot_data{k} = Data.data_u_dot{i}(:,1:p);
        v_dot_data{k} = Data.data_v_dot{i}(:,1:p);
        k = k + 1;
    end
end

% norm(u_dot_data{1}'*u_data{1} +(u_dot_data{1}'*u_data{1})','fro' )
% norm(u_dot_data{2}'*u_data{2} +(u_dot_data{2}'*u_data{2})','fro' )

% Prepare data for computing local coordinates
u_data_loc = cell(1,n);
v_data_loc = cell(1,n);
u_dot_data_loc = cell(1,n);
v_dot_data_loc = cell(1,n);

[U,Pu] = maxvol(u_data{1});
[V,Pv] = maxvol(v_data{1});
for i = 1:2
    u_data_loc{i} = Pu*u_data{i};
    v_data_loc{i} = Pv*v_data{i};
    u_dot_data_loc{i} = Pu*u_dot_data{i};
    v_dot_data_loc{i} = Pv*v_dot_data{i};
end


% Make plot
E_lag = [];
E_herm = [];
for i = 1:21
    t = points(i);
    I_lag = Interpolate_Gr(Ias, u_data,t, 'normal_lag');
    I_herm = Interpolate_Gr(Ias, u_data,t, 'normal_herm',u_dot_data);
    
    P = Data.data_u{i}(:,1:p)*Data.data_u{i}(:,1:p)';

    E_lag(i) = norm(P - I_lag*I_lag','fro');
    E_herm(i) = norm(P - I_herm*I_herm','fro');
end

f = figure
f.Position = [40,800,1200*5/6*1/2,650*5/6];
plot(points(1:21),E_lag)
hold on 
plot(points(1:21),E_herm)
title("Interpolation in normal coordinates")
legend("Lagrange","Hermite")
fontsize(f,15,"pixels")

% Make plot, local
E_lag = [];
E_herm = [];
for i = 1:21
    t = points(i);
    I_lag = Pu'*Interpolate_Gr(Ias, u_data_loc,t, 'local_lag');
    I_herm = Pu'*Interpolate_Gr(Ias, u_data_loc,t, 'local_herm',u_dot_data_loc);
    
    P = Data.data_u{i}(:,1:p)*Data.data_u{i}(:,1:p)';

    E_lag(i) = norm(P - I_lag*I_lag','fro');
    E_herm(i) = norm(P - I_herm*I_herm','fro');
end

f = figure
f.Position = [40,800,1200*5/6*1/2,650*5/6];
plot(points(1:21),E_lag)
hold on 
plot(points(1:21),E_herm)
title("Interpolation in local coordinates")
legend("Lagrange","Hermite")
fontsize(f,15,"pixels")





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