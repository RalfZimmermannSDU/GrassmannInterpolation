% Create snapshots
function FN_create_snapshots(points,num_time_pts)
    
    
    [~,n] = size(points);
    M = matrix_tools();
    
    nx = 1024;
    nt = 1001;
    
    
    Data = cell(1,n);
    Data_ph = cell(1,n);
    Data_dot = cell(1,n);
    
    data_u = cell(1,n);
    data_v = cell(1,n);
    
    data_u_dot = cell(1,n);
    data_v_dot = cell(1,n);
    
    for i = 1:n
        Y = FN_full_model(points(i),num_time_pts);
        %Y = Y / norm(Y);
        Data{i} = Y;
        
        %data_u{i} = Data{i}(1:nx,:);
        %data_v{i} = Data{i}(nx+1:end,:);
    
        h = 0.0000001;
        Ymh = FN_full_model(points(i)-h,num_time_pts);
        %Ymh = Y;
        Yph = FN_full_model((points(i)+h),num_time_pts);
        %Yph = Yph / norm(Yph);
        Data_dot{i} = (Yph - Ymh) / (2*h);
        % 
        Data_ph{i} = Yph;
        % % Compute u,v,u_dot and v_dot
        %[Up,Sp,~] = svd(Y(1:nx,1:nt),'econ');
        %[Vp,~,~] = svd(Y(nx+1:end,1:nt));
        %data_u{i} = Up(:,1:p);
        %Sp = Sp(1:p,1:p);
        %data_v{i} = Vp(:,1:p);
        % 

        %[Utemp,Stemp,~] = svd(Yph(1:nx,1:nt),'econ');
        %data_u_dot{i} = (Utemp(1:nx,1:p) - data_u{i}) / h;
        %data_u_dot{i} = M.dSVD(Y(1:nx,1:nt),p,Data_dot{i}(1:nx,1:nt));
        %data_v_dot{i} = M.dSVD(Data{i}(nx+1:end,1:nt),p,Data_dot{i}(nx+1:end,1:nt));
        disp("Computed data for Ia = " + num2str(points(i)))
    end
    % csvwrite('d1.csv',Data{1})
    % csvwrite('d2.csv',Data{2})
    % csvwrite('dd1.csv',Data_dot{1})
    % csvwrite('dd2.csv',Data_dot{2})

    %eval(['save snapshots_FN_model/snapshot_N_', num2str(n),'.mat Data Data_dot data_u data_v data_u_dot data_v_dot',' -v7.3']);   
    eval(['save snapshots_FN_model/snapshot_N_', num2str(n),'.mat Data Data_ph Data_dot',' -v7.3']);   
end
