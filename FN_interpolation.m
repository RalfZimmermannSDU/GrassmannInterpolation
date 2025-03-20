% FitzHugh-Nagumo interpolation experiment
clear 
close all

% Flags

createsnapshots = 0; % create new snapshots? 

visuals = 0; % show phase-space plots of the full system?
svd_plots = 0; % show svd and RIC plots for the snapshot data?
pmor = 1; % perform parametric model order reduction


Mat = matrix_tools(); % Import various Grassmann functions. 


%% Snapshots
% Create plots of the system
points = [0.03 0.04 0.05 0.06 0.07];


[~,n] = size(points);
Data = cell(1,n);
derivData = cell(1,n);
if createsnapshots
    for i = 1:n
        Data{i} = FN_full_model(points(i));
    end

    % Store data 
    %eval(['save snapshots_FN_model/','data','.mat Data']);'

    h = 0.0001;
    for i = 1:n
        Ymh = FN_full_model(points(i)-h);
        Y = Data{i};
        Yph = FN_full_model(points(i)+h);
        derivData{i} = (Yph - Ymh) / (2*h);
    end
    % Store derivative data 
    %eval(['save snapshots_FN_model/','derivative_data','.mat Data']);
end


Data = load('snapshots_FN_model/data_highres.mat');
d_Data = load('snapshots_FN_model/data_highres.mat');

Data = Data.Data; 
d_Data = d_Data.Data; 
% 
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
    p = 7;
    Y = Data{1};
    [nx,nt] = size(Y);
    
    Grassmann_points_u = cell(1,n);
    Grassmann_points_v = cell(1,n);
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
    
    % Friday: Test interpolated basis!
    % Ia \in [0.03, 0.07]
    Ia = 0.065;
    norm(Mat.ExpG(Grassmann_points_u{1},LagrangeInt(Ia,points,tangent_data_u)) - Grassmann_points_u{1})
    norm(Mat.ExpG(Grassmann_points_v{1},LagrangeInt(Ia,points,tangent_data_v)) - Grassmann_points_v{1})


    %Y = FN_reduced_model(Mat.ExpG(Grassmann_points_u{1},LagrangeInt(Ia,points,tangent_data_u)),Mat.ExpG(Grassmann_points_v{1},LagrangeInt(Ia,points,tangent_data_v)),Ia);
    
    % Compare full model with approximation
    %norm(FN_full_model(Ia) - Y,'fro')
    

    % Hermite interpolation
    % Ia  = 0.035
    
    Ia = 0.035;

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


    % Compute Delta_p1 and Delta_p2 via finite differences


    % Compute interpolant at Ia = 0.035
    
    [k,l] = size(Y);
    figure
    for j = 1:3:l
        plot(Y(j,:),Y(j + k/2,:),'r');
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