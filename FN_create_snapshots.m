% Create snapshots
clear 
%points = 0.03:0.01:0.12;
points = 0.03:0.001:0.12;
%points = 0.03;
[~,n] = size(points);
M = matrix_tools();

nx = 1024;
nt = 1001;
p = 12; % The user can truncate further later on

Data = cell(1,n);
Data_dot = cell(1,n);

data_u = cell(1,n);
data_v = cell(1,n);

data_u_dot = cell(1,n);
data_v_dot = cell(1,n);

for i = 1:n
    Data{i} = FN_full_model(points(i));

    h = 0.00001;
    Ymh = FN_full_model(points(i)-h);$
    Yph = FN_full_model(points(i)+h);
    Data_dot{i} = (Yph - Ymh) / (2*h);

    % Compute u,v,u_dot and v_dot
    [Up,~,~] = svd(Data{i}(1:nx,1:nt));
    [Vp,~,~] = svd(Data{i}(nx+1:end,1:nt));
    data_u{i} = Up(:,1:p);
    data_v{i} = Vp(:,1:p);

    data_u_dot{i} = M.dSVD(Data{i}(1:nx,1:nt),p,Data_dot{i}(1:nx,1:nt));
    data_v_dot{i} = M.dSVD(Data{i}(nx+1:end,1:nt),p,Data_dot{i}(nx+1:end,1:nt));
    disp("Computed data for Ia = " + num2str(points(i)))
end

eval(['save snapshots_FN_model/snapshot_N_', num2str(n),'.mat data_u data_v data_u_dot data_v_dot']);


