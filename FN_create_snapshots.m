% Create snapshots
function FN_create_snapshots(p,points,num_time_pts)
    
    
    [~,n] = size(points);
    M = matrix_tools();
    
    nx = 1024;
    nt = 1001;
    
    
    Data = cell(1,n);
    Data_dot = cell(1,n);
    
    data_u = cell(1,n);
    data_v = cell(1,n);
    
    data_u_dot = cell(1,n);
    data_v_dot = cell(1,n);
    
    for i = 1:n
        Data{i} = FN_full_model(points(i),num_time_pts);
        
        data_u{i} = Data{i}(1:nx,:);
        data_v{i} = Data{i}(nx+1:end,:);
    
        h = 0.00001;
        Ymh = FN_full_model(points(i)-h,num_time_pts);
        Yph = FN_full_model(points(i)+h,num_time_pts);
        Data_dot{i} = (Yph - Ymh) / (2*h);
        % 
        % % Compute u,v,u_dot and v_dot
        [Up,~,~] = svd(Data{i}(1:nx,1:nt));
        [Vp,~,~] = svd(Data{i}(nx+1:end,1:nt));
        data_u{i} = Up(:,1:p);
        data_v{i} = Vp(:,1:p);
        % 
        data_u_dot{i} = M.dSVD(Data{i}(1:nx,1:nt),p,Data_dot{i}(1:nx,1:nt));
        data_v_dot{i} = M.dSVD(Data{i}(nx+1:end,1:nt),p,Data_dot{i}(nx+1:end,1:nt));
        disp("Computed data for Ia = " + num2str(points(i)))
    end
    
    eval(['save snapshots_FN_model/snapshot_N_', num2str(n),'.mat data_u data_v data_u_dot data_v_dot',' -v7.3']);   

end
