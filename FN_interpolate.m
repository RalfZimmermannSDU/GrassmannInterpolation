% Interpolate the FN system
function FN_interpolate(p,Data,Data_ref,points,m,t0_glob,t1_glob,h,h2, maxsteps)

M = matrix_tools();

Intervals = t0_glob:h:t1_glob;
points = t0_glob:h2:t1_glob;

members = h/h2 + 1;

E_lag_loc = [];
E_herm_loc = [];
E_lag_norm = [];
E_herm_norm = [];

[~,mm] = size(points);

data_norms_l = [];
data_norms_r = [];

P_data_norms_l = [];
P_data_norms_r = [];

for j = 1:((t1_glob - t0_glob)/h)
    I = zeros(1,mm);
    t0 = Intervals(j);
    t1 = Intervals(j+1);

    Ias = [t0 t1];
    
    I = ismembertol(t0_glob:0.001:t1_glob,Ias,10e-10);
    n = sum(I);
    %[~,m] = size(points);
    
    % Extract data
    u_data = cell(1,n);
    v_data = cell(1,n);
    u_dot_data = cell(1,n);
    u_dot_data_hor = cell(1,n);
    v_dot_data_hor = cell(1,n);

    v_dot_data = cell(1,n);
    
    k = 1;
    for i = 1:m
        if I(i)
            u_data{k} = Data.data_u{i}(:,1:p);
            v_data{k} = Data.data_v{i}(:,1:p);
            u_dot_data{k} = Data.data_u_dot{i}(:,1:p);
            v_dot_data{k} = Data.data_v_dot{i}(:,1:p);

            % Compute the horizontal lifts
            u_dot_data_hor{k} = (u_dot_data{k}*u_data{k}'+u_data{k}*u_dot_data{k}')*u_data{k};
            v_dot_data_hor{k} = (v_dot_data{k}*v_data{k}'+v_data{k}*v_dot_data{k}')*v_data{k};
            
            % Check if the horizontal lifts are indeed horizontal
            norm(u_data{1}'*u_dot_data_hor{1},'fro')
            k = k + 1;

        end
        Data.data_u{i} = Data.data_u{i}(:,1:p);
    end
    
    true_data = cell(1,members);
    for k = 1:members
        true_data{k} = Data_ref.data_u{(members-1)*(j-1) + k}(:,1:p);
    end
    u_data_loc = cell(1,n);
    v_data_loc = cell(1,n);
    u_dot_data_loc = cell(1,n);
    v_dot_data_loc = cell(1,n);
    
    Pus = cell(1,2);
    Pvs = cell(1,2);

    for s = 1:2
        [~,Pus{s}] = maxvol(u_data{s},maxsteps);
        [~,Pvs{s}] = maxvol(v_data{s},maxsteps);
        for i = 1:2
            u_data_loc{i} = Pus{s}*u_data{i};
            v_data_loc{i} = Pvs{s}*v_data{i};
            %u_dot_data_loc{i} = Pu*u_dot_data{i};
            %v_dot_data_loc{i} = Pv*v_dot_data{i};
            conds_u(2*(s-1) + i) = cond(u_data_loc{i}(1:p,1:p),'fro');
            conds_v(2*(s-1) + i) = cond(v_data_loc{i}(1:p,1:p),'fro');
        end
    end
    
    for i = 1:2
        [~,iu] = min(max(conds_u(1:2),max(conds_u(3:4))));
        [~,iv] = min(max(conds_v(1:2),max(conds_v(3:4))));
    end
    
    Pu = Pus{iu};
    Pv = Pvs{iv};
    for i = 1:2
        u_data_loc{i} = Pu*u_data{i};
        v_data_loc{i} = Pv*v_data{i};
        u_dot_data_loc{i} = Pu*u_dot_data{i};
        v_dot_data_loc{i} = Pv*v_dot_data{i};
    end


    % Store before / after upper p x p block norm (Only data for the u - part is stored)
    data_norms_l(j) = norm(eye(p) / u_data{1}(1:p,1:p),'fro');
    data_norms_r(j) = norm(eye(p) / u_data{2}(1:p,1:p),'fro');
    P_data_norms_l(j) = norm(eye(p) / u_data_loc{1}(1:p,1:p),'fro');
    P_data_norms_r(j) = norm(eye(p) / u_data_loc{2}(1:p,1:p),'fro');

    [EL,EH] = error_on_interval_loc(u_data_loc,u_dot_data_loc,true_data,Pu,t0,t1,h2);
    E_lag_loc = [E_lag_loc(1:end-1) EL];
    E_herm_loc = [E_herm_loc(1:end-1) EH];
    
    [EL,EH] = error_on_interval(u_data,u_dot_data_hor,true_data,t0,t1,h2);
    E_lag_norm = [E_lag_norm(1:end-1) EL];
    E_herm_norm = [E_herm_norm(1:end-1) EH];
end

Left = data_norms_l';
Right = data_norms_r';
P_Left = P_data_norms_l';
P_Right = P_data_norms_r';

T = table(Left,P_Left,Right,P_Right);
disp(T);

[~,m] = size(E_lag_loc);
f = figure;
f.Position = [40,800,1200*5/6,650*5/6];
%points = 0.03:0.0001:0.12;
subplot(1,2,1)
plot(points(1:m),E_lag_loc)
hold on
plot(points(1:m),E_lag_norm,'--')
xlabel("I_a")
ylabel("Rel. error")
title("Error (Lagrange)")
legend("MV coords","Normal coords")



subplot(1,2,2)
plot(points(1:m),E_herm_loc)
hold on
plot(points(1:m),E_herm_norm,'--')
xlabel("I_a")
ylabel("Rel. error")
title("Error (Hermite)")
legend("MV coords","Normal coords")
fontsize(f,18,"pixels")
sgtitle("Relative interpolation errors")

exportgraphics(f,"experiment_2.png","Resolution",300);

end


function [E_lag,E_herm] = error_on_interval_loc(Data,deriv_data,true_data,Pd,t0,t1,h)
    ts = t0:h:t1;
    [~,m] = size(ts);
    Ias = [t0 t1];
    E_lag = [];
    E_herm = [];
    % feas_l = [];
    % feas_h = [];
    for i = 1:m
        t = ts(i);
        I_lag = Pd'*Interpolate_Gr(Ias, Data,t, 'local_lag');
        I_herm = Pd'*Interpolate_Gr(Ias, Data,t, 'local_herm',deriv_data);
        % feas_l(i) = norm(I_lag'*I_lag - eye(8),'fro');
        % feas_h(i) = norm(I_herm'*I_herm - eye(8),'fro');
        P = true_data{i}*true_data{i}';
        NP = norm(P,"fro");
    
        E_lag(i) = norm(P - I_lag*I_lag','fro')/NP;
        E_herm(i) = norm(P - I_herm*I_herm','fro')/NP;
    end
    
    % figure
    % plot(feas_h)
end

function [E_lag,E_herm] = error_on_interval(Data,deriv_data,true_data,t0,t1,h)
    ts = t0:h:t1;
    [~,m] = size(ts);
    Ias = [t0 t1];
    E_lag = [];
    E_herm = [];
    %feas_l = [];
    %feas_h = [];
    for i = 1:m
        t = ts(i);
        I_lag = Interpolate_Gr(Ias, Data,t, 'normal_lag');
        I_herm = Interpolate_Gr(Ias, Data,t, 'normal_herm',deriv_data);
        feas_l(i) = norm(I_lag'*I_lag - eye(8),'fro');
        feas_h(i) = norm(I_herm'*I_herm - eye(8),'fro');
        P = true_data{i}*true_data{i}';
        NP = norm(P,"fro");
        
        E_lag(i) = norm(P - I_lag*I_lag','fro')/NP;
        E_herm(i) = norm(P - I_herm*I_herm','fro')/NP;
    end
    %figure
    %plot(feas_l)
end