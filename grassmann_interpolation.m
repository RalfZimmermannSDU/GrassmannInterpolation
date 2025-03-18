 %% Grassmann interpolation

% Run interpolation experiments on the Grassmann manifold Gr(n,p):
%
% Curve 1: c(t) = B * A(t), A(t) = expm((t+1)^2*S), S \in skew(p).
% Curve 2: c(t) = B(t) * A(t), A(t) as above and B(t) = B0*(1-t) + B1*t.
% Curve 3: c(t) = B / A(t), A(t) as above.
% Curve 4: c(t) = B(t) / A(t), A(t) as above, B(t) = B * t^2 + 2 * B * t.
% Curve 5: c(t) = B * A(t), A(t) = expm((t+1)^2*S) / || expm((t+1)^2*S)||, 
%          S \in skew(p). Only for Lagrange
%
clear
rng(123123,'twister')
close all
% Genreal experiment parameters
lt0 = 0; lt1 = 10; % Time interval for Lagrange wave plot.
t0 = 0; t1 = 1; % Time interval for Hermite interpolation.
n = 50; p = 22; % Manifold dimension.
whatCurve = 1; % Which curve to use for the experiment.

% Interpolation paramters
lagrange = 1; % Run Lagrange interpolation experiment.
mm = 11; % Number of Lagrange interpolation points on [lt0,lt1].
MM = 501; % Number of Lagrange interpolants on [lt0,lt1].

useChebyshev = 0; % Use Chebyshev  points instead of equidistant points.
hermite = 1; % Run Hermite interpolation experiment (also showing Lagrange).
m = 51; % Number of interpolations to perform in [t0,t1].

M = matrix_tools(); % Import various Grassmann functions. 

% Curves
if whatCurve == 1
    % Create data from the curve t -> B * exp((t+1)^2*S), for S \in
    % \skew(p).
    S = rand(p);
    S = S / norm(S,'fro');
    S = 0.5 * (S - S');
    B = rand(n-p,p);
    B = B / norm(B,'fro');
    
    c  = @(t) B*expm((t+1)^2*S);
    dc = @(t) 2*B*(t+1)*S*expm((t+1)^2*S);
elseif whatCurve == 2
    % Create data from the curve t -> B(t) * exp((t+1)^2*S), for S \in
    % \skew(p).
    % and B(t) is a line between two points in R^((n-p) x p).
    S = rand(p);
    S = 0.5 * (S - S');
    S = S / norm(S,'fro');
    B0 = rand(n-p,p);
    B0 = B0 / norm(B0,'fro');
    B1 = rand(n-p,p);
    B1 = B1 / norm(B1,'fro');
    
    b = @(t) B0*(1-t) + t*B1;
    c = @(t) b(t) * expm((t+1)^2*S);
    dc = @(t) (B1 - B0) * expm((t+1)^2*S) + 2*b(t)*(t+1)*S*expm((t+1)^2*S);
elseif whatCurve == 3
    % Create data from the curve t -> B / exp((t+1)^2*S), for S \in
    % \skew(p).
    % Consider the inverse case of curve 1.
    S = rand(p);
    S = 0.5 * (S - S');
    
    B = rand(n-p,p);
    B = B / norm(B,'fro');

    c  = @(t) B / expm((t+1)^2*S);
    dc = @(t) -2*B / expm((t+1)^2*S) *(t+1)*S*expm((t+1)^2*S)/expm((t+1)^2*S);

elseif whatCurve == 4
    % Create data from the curve t -> B(t) / exp((t+1)^2*S), for S \in
    % \skew(p).
    % Let B(t) be a polynomial.
    S = rand(p);
    S = 0.5 * (S - S');
    B = rand(n-p,p);
    B = B / norm(B,'fro');

    b = @(t) B * t^2 + 2*B*t;

    c  = @(t) b(t) / expm((t+1)^2*S);
    dc = @(t) (2*B*t + 2*B)/ expm((t+1)^2*S) - 2*b(t) / expm((t+1)^2*S) *(t+1)*S*expm((t+1)^2*S)/expm((t+1)^2*S);
elseif whatCurve == 5
    % Create data from the curve  
    % t -> B / A(t), A(t) = expm((t+1)^2*S) / || expm((t+1)^2*S)||, 
    %      with S \in skew(p). Only for Lagrange
    S = rand(p);
    S = S / norm(S,'fro');
    S = 0.5 * (S - S');
    B = rand(n-p,p);
    B = B / norm(B,'fro');

    c  = @(t) B / (expm((t+1)^2*S)/norm(expm((t+1)^2*S),'fro'));
    dc = @(t) 0; % Not relevant for Lagrange.
else
    error("No such curve exist!")
end

% Check that we computed the derivative dc correctly using FD:
t = 0.0000001;
g = (c(t+t0) - c(t0)) / t;
h = dc(t0);
acct0 = norm(g-h,'fro');

g = (c(t+t1) - c(t1)) / t;
h = dc(t1);
acct1 = norm(g-h,'fro');

disp("The derivative at t0 has FD error " + num2str(acct0))
disp("The derivative at t1 has FD error " + num2str(acct1))

% The data for Lagrange interpolation is structured slightly different 
% than the data for Hermite interpolation. 

% The data is generated in local coordinartes, then mapped to Gr(n,p).
if lagrange
    if useChebyshev
        % Genreate Chebyshev nodes.
    else
        % Equidistant nodes.
        ts = linspace(lt0,lt1,mm);
    end
    
    % Generate and map 
    local_ys = cell(1,mm);
    ys = cell(1,mm);
    extsum = 0; % To check data validity.
    for i = 1:mm
        local_ys{i} = c(ts(i)); % Data in local cooridnates
        ys{i} = M.ParamG(local_ys{i}); % Data on Grassmann
        extsum = extsum + norm(ys{i}^2 - ys{i},'fro');
    end
    extsum = extsum / mm; % Average feasibility. Should be small. 
    
    % To compare to interpolation in Riemannian coordinates, generate
    % tangent space data by generating the Stiefel representation of each
    % data point. Choose a tangent space (the one for the leftmost point
    % by default) and map the data. 
    tangent_data = cell(1,mm);
    [U,~,~] = svd(ys{1});
    U = U(:,1:p); % Stiefel representation
    for i = 1:mm
        [V,~,~] = svd(ys{i});
        V = V(:,1:p); % Stiefel representation of Grassmann data
        tangent_data{i} = M.LogG(U,V); % Data is mapped to tangent space
    end
end

% Lagrange wave plot
if lagrange
    
    int_err_tan = []; 
    int_err_loc = [];
    feas_tan = [];
    feas_loc = [];

    tsInt = linspace(lt0,lt1,MM);
    % Lagrange in normal coordinates (in the tangent space)
    for i = 1:MM
        Int = LagrangeInt(tsInt(i),ts,tangent_data);
        Int_S_rep = M.ExpG(U,Int); % U is the anchor, int is the interpolant
        Int_Proj = Int_S_rep*Int_S_rep';
        int_err_tan(i) = norm(Int_Proj - M.ParamG(c(tsInt(i))),'fro'); % Interpolation error
        feas_tan(i) = norm(Int_Proj*Int_Proj - Int_Proj,'fro'); % Feasibility
    end
    
    % Lagrange in local coordinates
    for i = 1:MM
        Int = LagrangeInt(tsInt(i),ts,local_ys);
        Int_Proj = M.ParamG(Int);
        int_err_loc(i) = norm(Int_Proj - M.ParamG(c(tsInt(i))),'fro');
        feas_loc(i) = norm(Int_Proj*Int_Proj - Int_Proj,'fro'); % Feasibility
    end

    f = figure; 
    f.Position = [40,800,1200*5/6,650*5/6];
    subplot(1,2,1)
    plot(tsInt,int_err_tan,'LineWidth',2)
    hold on
    plot(tsInt,int_err_loc,'LineWidth',2)
    legend('Lagrange normal','Lagrange local')
    title("Interpolation error")
    
    subplot(1,2,2)
    semilogy(tsInt,feas_tan,'LineWidth',2)
    hold on
    semilogy(tsInt,feas_loc,'LineWidth',2)
    
    legend('Lagrange normal','Lagrange Local')
    title("Feasibility ||P^2 - P||_F")
    fontsize(f,15,"pixels")
    exportgraphics(f,'output_lagrange_wave.png','Resolution',300)
end


% Hermite one-wave plot 
if hermite
    % Generate data at the start- and endpoint
    % All data is in this case given from the projector perspective
    x0 = M.ParamG(c(t0));
    x1 = M.ParamG(c(t1));
    v0 = M.dParamG(c(t0),dc(t0));
    v1 = M.dParamG(c(t1),dc(t1));
    
    % % Check that we recover c(tt) and dc(tt)
    % norm(c(t1) - M.LocalCoordG(x1,n,p))
    % norm(dc(t1) - M.dLocalCoordG(x1,v1,n,p))

    % Interpolate on [t0,t1]
    m = 51;
    t = linspace(t0,t1,m);

    % Local coordinate data. 
    Locx0 = M.LocalCoordG(x0,n,p);
    Locx1 = M.LocalCoordG(x1,n,p);
    dLocx0 = M.dLocalCoordG(x0,v0,n,p);
    dLocx1 = M.dLocalCoordG(x1,v1,n,p);

    % Tangent space data. Take x1 as anchor.
    % Construct n x p curve for QR decomposition-based Stiefel parameteization
    cprime = @(t) [eye(p);c(t)];
    [q0,~] = qr(cprime(t0),'econ');
    [q1,~] = qr(cprime(t1),'econ');
    
    % % Sanity check: Does the projector x0 match q0 * q0'?
    % norm(x0 - q0*q0')
    % norm(x1 - q1*q1')

    xi = M.LogG(q1,q0); % With x1 as base, compute xi
    % % Sanity check: Is xi the correct Riemannian logarithm?
    % norm(q0*q0' - M.ExpG(q1,xi,1)*M.ExpG(q1,xi,1)','fro')
    
    % The tangent vectors are lifted to Stiefel tangents
    hor_v0 = v0*q0;
    hor_v1 = v1*q1;
    % % Sanity check: Are the lifts indeed true Stiefel tangents?
    % norm(hor_v0'*q0 + q0'*hor_v0,'fro')
    % norm(hor_v1'*q1 + q1'*hor_v1,'fro')
    h = 1e-4;
    DELTA_p = FDapprox(M, q0, q1, hor_v0, h); % This is a tangent vector at q1
    
    % % Sanity check: Is the FD approx a Stiefel tangent at q1?
    % norm(DELTA_p'*q1 + q1'*DELTA_p,'fro')
    
    a0 = @(t) 1 - 1/(t1-t0)^2*(t-t0)^2 + 2/(t1 - t0)^3*(t-t0)^2*(t-t1);
    a1 = @(t) 1/(t1-t0)^2 * (t-t0)^2 - 2/(t1-t0)^3*(t-t0)^2*(t-t1);
    b0 = @(t) (t-t0) - 1/(t1-t0)*(t-t0)^2 + 1/(t1-t0)*(t-t0)^2*(t-t1);
    b1 = @(t) 1/(t1-t0)^2*(t-t0)^2*(t-t1);

    % Tangent space curve
    mu = @(t) a0(t) * xi + b0(t) * DELTA_p + b1(t) * hor_v1;
    
    % Compare to Lagrange interpolation in local coordinates.
    ts_lag = [t0,t1];
    celly = cell(1,2);
    celly{1} = Locx0; celly{2} = Locx1;

    Ielag = [];
    Ieherm = [];
    Ieherm_tan = [];

    feas_loc_lag = [];
    feas_loc_herm = [];
    feas_norm_herm = [];
    
    % Compute the interpolants
    for i = 1:m
        ythermite = HermiteInterpol(Locx0,Locx1,dLocx0,dLocx1,t0,t1,t(i));
        ytlagrange = LagrangeInt(t(i),ts_lag,celly);
        Ieherm(i) = norm(M.ParamG( c(t(i)) ) - M.ParamG( ythermite ) ,'fro');
        Ielag(i) = norm(M.ParamG( c(t(i)) ) - M.ParamG( ytlagrange ) ,'fro');
        
        y_herm_tan = M.ExpG(q1,mu(t(i)))*M.ExpG(q1,mu(t(i)))';
        Ieherm_tan(i) = norm(y_herm_tan*y_herm_tan' - M.ParamG( c(t(i)) ),'fro');
        
        feas_loc_lag(i) = norm((M.ParamG(ytlagrange)*M.ParamG(ytlagrange) - M.ParamG(ytlagrange)),'fro');
        feas_loc_herm(i) = norm((M.ParamG(ythermite)*M.ParamG(ythermite) - M.ParamG(ythermite)),'fro');
        feas_norm_herm(i) = norm(y_herm_tan*y_herm_tan - y_herm_tan,'fro');
    end
    

    f = figure; 
    f.Position = [40,800,1200*5/6,650*5/6];
    subplot(1,2,1)
    plot(t,Ieherm,'LineWidth',2)
    hold on
    plot(t,Ieherm_tan,'LineWidth',2)
    plot(t,Ielag,'LineWidth',2)
    legend('Hermite local','Hermite normal','Lagrange local')
    title("Interpolation error")
    
    subplot(1,2,2)
    semilogy(t,feas_loc_herm,'LineWidth',2)
    hold on
    semilogy(t,feas_norm_herm,'LineWidth',2)
    semilogy(t,feas_loc_lag,'LineWidth',2)

    legend('Hermite local','Hermite normal','Lagrange local')   
    title("Feasibility ||P^2 - P||_F")
    fontsize(f,15,"pixels")
    exportgraphics(f,'output_hermite_compare.png','Resolution',300)

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
function deriv = FDapprox(Man, y0,y1,dir,h)
    % Input:
    %   Man:   Manifold object
    %   y0,y1: Points 
    %   dir:   Direction in which the differnetial should be computed
    %   h:     Stepsize

    M = Man;
    deriv = ( M.LogG(y1,M.ExpG(y0,h*dir)) - M.LogG(y1,M.ExpG(y0,-h*dir)) ) / (2*h);
end