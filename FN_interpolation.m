% FitzHugh-Nagumo interpolation experiment
clear 
close all

% Flags

createsnapshots = 0; % create new snapshots? 

visuals = 0; % show phase-space plots of the full system?
svd_plots = 0; % show svd and RIC plots for the snapshot data?
pmor = 1; % perform parametric model order reduction

interpol_bases_plot = 0; % comare true and interpolated reduce bases 

Mat = matrix_tools(); % Import various Grassmann functions. 


%% Snapshots
% Create plots of the system
points = [0.03 0.04 0.05 0.06 0.07];
%points = [0.07 0.09 0.11];

[~,n] = size(points);
Data = cell(1,n);
derivData = cell(1,n);
if createsnapshots
    for i = 1:n
        Data{i} = FN_full_model(points(i));
    end

    % Store data 
    eval(['save snapshots_FN_model/','data_shift','.mat Data']);

    h = 0.0001;
    for i = 1:n
        Ymh = FN_full_model(points(i)-h);
        Y = Data{i};
        Yph = FN_full_model(points(i)+h);
        derivData{i} = (Yph - Ymh) / (2*h);
    end
    % Store derivative data 
    eval(['save snapshots_FN_model/','derivative_data_shift','.mat Data']);
end


% Data = load('snapshots_FN_model/data_highres.mat');
% d_Data = load('snapshots_FN_model/derivative_data.mat');
% 
% Data = Data.Data; 
% d_Data = d_Data.Data; 

Data = load('snapshots_FN_model/data_highres.mat');
d_Data = load('snapshots_FN_model/derivative_data.mat');
Data = Data.Data; 
d_Data = d_Data.Data; 

% for i = 1:n
%     rank(Data{i})
% end
% 

%% SVD and RIC plots
% Create plots showing the decay of the singular values and relative
% information content plots.
if svd_plots
    f = figure;
    for i = 1:n
        S = svd(Data{i});
        
        semilogy(S(1:100)/sum(S),'x');
        hold on
    end
    fontsize(f,15,"pixels")
    legend("I_a = " + num2str(points(1)),"I_a = " + num2str(points(2)),"I_a = " + num2str(points(3)),...,
        "I_a = " + num2str(points(4)),"I_a = " + num2str(points(5)))
    xlabel("Position of singular value in \Sigma")
    title("Singular values full snapshot matrices")
    
    % RIC curve
    M = 20;
    f = figure;
    for i = 1:5
        information = zeros(1,M);
        S = svd(Data{i});
        for j = 2:M+1
            information(j) = sum(S(1:j-1)) / sum(S);
        end
        plot(0:M,information,'LineWidth',2)
        hold on
    end
    
    fontsize(f,15,"pixels")
    legend("I_a = " + num2str(points(1)),"I_a = " + num2str(points(2)),"I_a = " + num2str(points(3)),...,
        "I_a = " + num2str(points(4)),"I_a = " + num2str(points(5)))
    xlabel("Number of singular values")
    ylabel("RIC")
    title("RIC curves for the full snapshot matrices")
    
    % Select u and v component
    % u component 
    f = figure;
    f.Position = [40,800,1200*5/6,650*5/6];
    D = Data{1};
    [nx,nt] = size(D);
    subplot(1,2,1)
    for i = 1:n
        D = Data{i};
        S = svd(D(1:nx/2,1:nt));
        
        semilogy(S(1:100)/sum(S),'x');
        hold on
    end
    fontsize(f,15,"pixels")
    legend("I_a = " + num2str(points(1)),"I_a = " + num2str(points(2)),"I_a = " + num2str(points(3)),...,
        "I_a = " + num2str(points(4)),"I_a = " + num2str(points(5)))
    xlabel("Position of singular value in \Sigma")
    title("u-component")
    
    % v component 
    %f = figure;
    subplot(1,2,2)
    for i = 1:n
        D = Data{i};
        S = svd(D(nx/2+1:end,1:nt));
        
        semilogy(S(1:100)/sum(S),'x');
        hold on
    end
    fontsize(f,15,"pixels")
    legend("I_a = " + num2str(points(1)),"I_a = " + num2str(points(2)),"I_a = " + num2str(points(3)),...,
        "I_a = " + num2str(points(4)),"I_a = " + num2str(points(5)))
    xlabel("Position of singular value in \Sigma")
    title("v-component")
    
    sgtitle("Component-wise singular values")
    
    f = figure
    f.Position = [40,800,1200*5/6,650*5/6];
    
    subplot(1,2,1)
    % RIC curves
    % u-component
    M = 20;
    for i = 1:5
        information = zeros(1,M);
        Y = Data{i};
        S = svd(Y(1:nx/2,:));
        for j = 2:M+1
            information(j) = sum(S(1:j-1)) / sum(S);
        end
        plot(0:M,information,'LineWidth',2)
        hold on
    end
    fontsize(f,15,"pixels")
    legend("I_a = " + num2str(points(1)),"I_a = " + num2str(points(2)),"I_a = " + num2str(points(3)),...,
        "I_a = " + num2str(points(4)),"I_a = " + num2str(points(5)))
    xlabel("Number of singular values")
    ylabel("RIC")
    title("u-component")
    
    subplot(1,2,2)
    % v-component
    M = 20;
    for i = 1:5
        information = zeros(1,M);
        Y = Data{i};
        S = svd(Y(nx/2+1:end,:));
        for j = 2:M+1
            information(j) = sum(S(1:j-1)) / sum(S);
        end
        plot(0:M,information,'LineWidth',2)
        hold on
    end
    fontsize(f,15,"pixels")
    legend("I_a = " + num2str(points(1)),"I_a = " + num2str(points(2)),"I_a = " + num2str(points(3)),...,
        "I_a = " + num2str(points(4)),"I_a = " + num2str(points(5)))
    xlabel("Number of singular values")
    ylabel("RIC")
    title("v-component")
    
    sgtitle("RIC curves for the full snapshot matrices")
end

%% Interpolation 

% Perform parametric model order reduction of the FN system

if pmor
    % Generate the points on each of the two Grassmann manifolds 
    p = 8;
    Y = Data{1};
    [nx,nt] = size(Y);
    
    Grassmann_points_u = cell(1,n);
    Grassmann_points_v = cell(1,n);
    
    % Use fewer snapshots
    %points = [0.03, 0.05, 0.07];

    for i = 1:n
        Y = Data{i};
    
        [Uup,~,~] = svd(Y(1:nx/2,1:nt));
        [Uvp,~,~] = svd(Y(nx/2+1:end,1:nt));
        Grassmann_points_u{i} = Uup(:,1:p);
        Grassmann_points_v{i} = Uvp(:,1:p);
    end
    % 
    % for i = 1:n
    %     Uu = Grassmann_points_u{i};
    %     Uu = Uu(1:p,1:p);
    %     cond(Uu)
    % end
    
    % Interpolate in the tangent space via Lagrange
    tangent_data_u = cell(1,n);
    tangent_data_v = cell(1,n);
    U = Grassmann_points_u{1};
    W = Grassmann_points_v{1};
    for i = 1:n
        V = Grassmann_points_u{i};
        Z = Grassmann_points_v{i};
        tangent_data_u{i} = Mat.LogG(U,V); % Data is mapped to tangent space
        tangent_data_v{i} = Mat.LogG(W,Z); % Data is mapped to tangent space
    end
    
    % Test interpolated basis!
    % Ia \in [0.03, 0.07]
    Ia = 0.035;
    Int1 = Interpolate_Gr(points, Grassmann_points_v,Ia, 'normal_lag');
    norm(Mat.ExpG(Grassmann_points_u{1},LagrangeInt(Ia,points,tangent_data_u)) - Grassmann_points_u{1})
    norm(Mat.ExpG(Grassmann_points_v{1},LagrangeInt(Ia,points,tangent_data_v)) - Grassmann_points_v{1})

    %YL = FN_reduced_model(Mat.ExpG(Grassmann_points_u{1},LagrangeInt(Ia,points,tangent_data_u)),Mat.ExpG(Grassmann_points_v{1},LagrangeInt(Ia,points,tangent_data_v)),Ia);
    
    % Compare full model with approximation
    %Y_full = FN_full_model(Ia);
    %e1 = norm(Y_full - YL,'fro');
    

    % Hermite interpolation
    % Ia  = 0.035
    
    % W1 = Data{2};
    % W2 = Data{3};
    % dW1 = d_Data{2};
    % dW2 = d_Data{3};
    % Data = cell(1,2);
    % d_Data = cell(1,2);
    % 
    % Data{1} = W1; Data{2} = W2;
    % d_Data{1} = dW1; d_Data{2} = dW2;

    % Grassmann data
    % Data{1} corresponds to  Ia = 0.03, Data{2} corresponds to Ia = 0.04 
    [U0,S0,V0] = svd(Data{1}(1:nx/2,:));
    [U1,S1,V1] = svd(Data{2}(1:nx/2,:));

    U1x0 = U0(:,1:p); % on Gr(n,p)
    U1x1 = U1(:,1:p); % on Gr(n,p)
    
    [U0,S0,V0] = svd(Data{1}(nx/2+1:end,:));
    [U1,S1,V1] = svd(Data{2}(nx/2+1:end,:));

    U2x0 = U0(:,1:p); % on Gr(n,p)
    U2x1 = U1(:,1:p); % on Gr(n,p)

    % | U1    |
    % |    U2 | is the reduced basis 

    % Derivative data
    %
    % Obtain derivatives for the four U factors
    [U1x0dot, Sdot, Vdot] = Mat.dSVD(Data{1}(1:nx/2,:),p,d_Data{1}(1:nx/2,:));
    [U1x1dot, Sdot, Vdot] = Mat.dSVD(Data{2}(1:nx/2,:),p,d_Data{2}(1:nx/2,:));

    [U2x0dot, Sdot, Vdot] = Mat.dSVD(Data{1}(nx/2+1:end,:),p,d_Data{1}(nx/2+1:end,:));
    [U2x1dot, Sdot, Vdot] = Mat.dSVD(Data{2}(nx/2+1:end,:),p,d_Data{2}(nx/2+1:end,:));
    
    % Asses quality 
    norm(U1x0dot'*U1x0+U1x0'*U1x0dot)
    norm(U1x1dot'*U1x1+U1x1'*U1x1dot)

    norm(U2x0dot'*U2x0+U1x0'*U2x0dot)
    norm(U2x1dot'*U2x1+U1x1'*U2x1dot)
    
    % Compute xi1 and xi2
    xi1 = Mat.LogG(U1x1,U1x0);
    xi2 = Mat.LogG(U2x1,U2x0);

    % Compute Delta_p1 and Delta_p2 via finite differences
    h = 0.0001;
    Delta_p1 = FDapprox(Mat, U1x0, U1x1, U1x0dot,h);
    Delta_p2 = FDapprox(Mat, U2x0, U2x1, U2x0dot, h);

    norm(Delta_p1'*U1x1+U1x1'*Delta_p1)
    norm(Delta_p2'*U2x1+U2x1'*Delta_p2)

    % Compute interpolant at Ia = 0.035
    D1 = HermiteInterpol(xi1,0,Delta_p1,U1x1dot,points(1),points(2),Ia);
    Q1 = Mat.ExpG(U1x1,D1);

    D2 = HermiteInterpol(xi2,0,Delta_p2,U2x1dot,points(1),points(2),Ia);
    Q2 = Mat.ExpG(U2x1,D2);
    
    derivs = cell(1,2);
    derivs{1} = U2x0dot;
    derivs{2} = U2x1dot;

    Int2 = Interpolate_Gr(points, Grassmann_points_v,Ia, 'normal_lag', derivs);

    %YH = FN_reduced_model(Q1,Q2,Ia);
    
    %e2 = norm(Y_full - YH,'fro');
    
    %
    % Interpolation in local coordinates
    %
    % Lagrange:

    % Preprocessing by maximum volume method 
    P_data_u = cell(1,n);

    [P_Gr_u1,Pu] = maxvol(Grassmann_points_u{1});
    
    mu = 0;
    for i = 1:n
        w =Pu*Grassmann_points_u{i};
        P_data_u{i} = w;

        mu = mu + cond(w(1:p,:));
    end
    disp("Average upper block condition number (u) "+ num2str(mu/n))

    P_data_v = cell(1,n);

    [P_Gr_v1,Pv] = maxvol(Grassmann_points_v{1});
    
    mu = 0;
    for i = 1:n
        w =Pv*Grassmann_points_v{i};
        P_data_v{i} = w;

        mu = mu + cond(w(1:p,:));
    end
    disp("Average upper block condition number (v) "+ num2str(mu/n))

    Q1 = Pu'*Interpolate_Gr(points, P_data_u,Ia, 'local_lag');
    Q2 = Pv'*Interpolate_Gr(points, P_data_v,Ia, 'local_lag');

    YL = FN_reduced_model(Q1,Q2,Ia);
    Y_full = FN_full_model(Ia);
    e3 = norm(Y_full - YL,'fro');


    % Hermite interpolation
    % Use the processed data from the Lagrange experiment
    
    % Compute the derivatives in local cooridnates 
    P_derivs_u = cell(1,2);
    P_derivs_u{1} = Pu*U1x0dot;
    P_derivs_u{2} = Pu*U1x1dot;

    P_derivs_v = cell(1,2);
    P_derivs_v{1} = Pu*U2x0dot;
    P_derivs_v{2} = Pu*U2x1dot;
    
    % % Map to local coordinates
    % for i = 1:2
    %     P_loc_derivs_u{i} = Mat.dLocalCoordG(P_data_u{i}, P_loc_derivs_u{i}, nx/2, p);
    %     P_loc_derivs_v{i} = Mat.dLocalCoordG(P_data_v{i}, P_loc_derivs_v{i}, nx/2, p);
    % end

    Q1 = Pu'*Interpolate_Gr(points(1:2), P_data_u,Ia, 'local_herm', P_derivs_u);
    Q2 = Pv'*Interpolate_Gr(points(1:2), P_data_v,Ia, 'local_herm', P_derivs_v);
    
    YH = FN_reduced_model(Q1,Q2,Ia);
    e4 = norm(Y_full - YH,'fro');
    
    [k,l] = size(YH);
    figure
    for j = 1:3:l
        plot(YH(j,:),YH(j + k/2,:),'r');
        hold on
    end
    xlabel('u(x,t)')
    ylabel('v(x,t)')
    title("Phase space, I_a = " + num2str(Ia))
end
%% Visuals
% Create plots of the system
if visuals
    for i = 1:n
        Y = Data{i};
        [k,l] = size(Y);
        figure
        for j = 1:3:l
            plot(Y(j,:),Y(j + k/2,:),'r');
            hold on
        end
        xlabel('u(x,t)')
        ylabel('v(x,t)')
        title("Phase space, I_a = " + num2str(points(i)))
    end
    
    f = figure;
    f.Position = [40,800,1200*5/6,650*5/6];
    subplot(1,2,1)
    %figure
    i = 1;
    Y = Data{i};
    [k,l] = size(Y);
    p1 = plot(Y(1,:),Y(1 + k/2,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',2);
    hold on
    for j = 4:3:l
        plot(Y(j,:),Y(j + k/2,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',2);
    end
    
    i = 3;
    Y = Data{i};
    [k,l] = size(Y);
    p2 = plot(Y(1,:),Y(1 + k/2,:),'Color',[0.9290 0.6940 0.1250],'LineWidth',2);
    hold on
    for j = 4:3:l
        plot(Y(j,:),Y(j + k/2,:),'Color',[0.9290 0.6940 0.1250],'LineWidth',2	);
    end
    
    i = 5;
    Y = Data{i};
    [k,l] = size(Y);
    p3 = plot(Y(1,:),Y(1 + k/2,:),'Color',[0.4940 0.1840 0.5560],'LineWidth',2);
    hold on
    for j = 4:3:l
        plot(Y(j,:),Y(j + k/2,:),'Color',[0.4940 0.1840 0.5560],'LineWidth',2);
    end
    
    xlabel('u(x,t)')
    ylabel('v(x,t)')
    title("Phase space, I_a = " + num2str(points(1))+ ', ' + num2str(points(3))+ ', '+ num2str(points(5))+ ', ' )
    legend([p1,p2,p3],{"I_a = "+num2str(points(1)),"I_a = "+num2str(points(3)),"I_a = "+num2str(points(5))})
    fontsize(f,15,"pixels")

    %exportgraphics(f,'FN_full_three_phase.png','Resolution',300)


    %f = figure;
    subplot(1,2,2)

    i = 1;
    Y = Data{i};
    
    [nx,nt] = size(Y);
    plot3(0*ones(1,nt),Y(1,:),Y(nx/2+1,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
    hold on
    x = 1/2^3:1/2^3:1;
    for j = 1:8
        plot3(x(j)*ones(1,nt),Y(128*j,:),Y(nx/2+128*j,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',2);
    end
    xlabel('x')
    ylabel('u(x,t)')
    zlabel('v(x,t)')
    title("Phase space, I_a = " + num2str(points(1)))

    fontsize(f,15,"pixels")
    exportgraphics(f,'FN_full.png','Resolution',300)
end

function interpol = LagrangeInt(x,xs,ys)
    n = length(xs);
    L = ones(1,n);
    interpol = 0;
    for j = 1:n
        for i = 1:n
            if i ~= j
                L(j) = L(j)*(x - xs(i)) / (xs(j)-xs(i));
            end
        end
    end
    for i = 1:n
        interpol = interpol + ys{i}*L(i);
    end
end

function y = HermiteInterpol(x0,x1,dx0,dx1,t0,t1,t)
    % Computes the Hermite Interpolant at t 
    % in a vector space.
    %   
    %   Input: Data x0, x1 in vector space format (ex. local coordinates)
    %          Derivative data dx0, dx1 in the same vector space
    %
    %   Output: Interpolated point at time t \in [t0,t1]
    %

    a0 = @(t) 1 - 1/(t1-t0)^2*(t-t0)^2 + 2/(t1 - t0)^3*(t-t0)^2*(t-t1);
    a1 = @(t) 1/(t1-t0)^2 * (t-t0)^2 - 2/(t1-t0)^3*(t-t0)^2*(t-t1);
    b0 = @(t) (t-t0) - 1/(t1-t0)*(t-t0)^2 + 1/(t1-t0)*(t-t0)^2*(t-t1);
    b1 = @(t) 1/(t1-t0)^2*(t-t0)^2*(t-t1);

    gamma = @(t) a0(t)*x0 + a1(t)*x1 + b0(t)*dx0 + b1(t) * dx1;

    y = gamma(t);
end
function deriv = FDapprox(Man, y0,y1,dir,h)
    % Input:
    %   Man:   Manifold object
    %   y0,y1: Points 
    %   dir:   Direction in which the differnetial should be computed
    %   h:     Stepsize

    M = Man;
    deriv = ( M.LogG(y1,M.ExpG(y0,h*dir)) - M.LogG(y1,M.ExpG(y0,-h*dir)) ) / (2*h);
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
