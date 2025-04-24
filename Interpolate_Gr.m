function Y = Interpolate_Gr(time_data, Data, t, routine_flag, Deriv_data)
% This function interpolates data on the Grassmann manifold sterming from a
% function f:[a,b] -> Gr(n,p)
%
% Input:
%       time_data:      Points in time for Data and Deriv_data as a list
%       Data:           Point data in Stiefel form on Gr(n,p) as a cell
%       Deriv_data:     Optional: Derivative data in Stiefel form as a cell
%       t:              Time where the interpolant is computed
%       routine_flag:   Can be 'normal_lag', 'normal_herm', 'local_lag' or
%                       'local_herm'
%       
% Output:
%       Y:              Interpolant
%

Mat = matrix_tools();

if nargin < 5 && isempty('Deriv_data')
    Deriv_data = cell(0);
end

if routine_flag == "normal_lag"
    % Preprocessing: Map data to the tangent space at the first data point
    [~,n] = size(Data);

    tangent_data = cell(1,n);
    
    U = Data{1}; % Base point
    for i = 1:n
        tangent_data{i} = Mat.LogG(U, Data{i});
    end

    % Interpolate
    Y = LagrangeInt(t, time_data, tangent_data);

    % Map back to manifold 
    Y = Mat.ExpG(U,Y);
end

if routine_flag == "normal_herm"
    % Preprocessing: Map data to the tangent space and compute quantities
    % use the rightmost point as reference.
    
    p = Data{1};
    q = Data{2};
    pdot = Deriv_data{1};
    qdot = Deriv_data{2};

    xi = Mat.LogG(q,p);
    
    h = 0.001;
    Delta_p = FDapprox(Mat, p, q, pdot, h);
    
    % Interpolate
    Y = HermiteInterpol(xi, 0, Delta_p, qdot, time_data(1), time_data(2), t);

    % Map back to manifold
    Y = Mat.ExpG(q,Y);

end

if routine_flag == "local_lag"
    % Preprocessing: Map to local cooridnates 
    [~,n] = size(Data);

    local_data = cell(1,n);
    
    U = Data{1}; 
    [m,p] = size(U);
    for i = 1:n
        % psi(U) = U2 / U1;
        local_data{i} = Mat.LocalCoordG(Data{i},m,p);
    end

    % Interpolate 
    Y = LagrangeInt(t, time_data, local_data);

    
    % Map back to manifold
    [Y,~] = qr([eye(p); Y],'econ');

    % R = chol(eye(p) + Y'*Y); % Note the transpose for consistency
    % Y = [eye(p); Y] / inv(R);
end

if routine_flag == "local_herm"
    % Preprocessing: Map to local coordinates
    p = Data{1};
    q = Data{2};
    pdot = Deriv_data{1};
    qdot = Deriv_data{2};

    [n,k] = size(p);
    
    lp = Mat.LocalCoordG(p,n,k);
    lq = Mat.LocalCoordG(q,n,k);
    
    lpdot = Mat.dLocalCoordG(p,pdot,n,k);
    lqdot = Mat.dLocalCoordG(q,qdot,n,k);
    
    % Interpolate
    Y = HermiteInterpol(lp, lq, lpdot, lqdot,time_data(1), time_data(2), t);
    
    % Map back to manifold
    R = chol(eye(k) + Y'*Y); % Note the transpose for consistency
    Y = [eye(k); Y] / R;
end

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
function y = HermiteInterpol(x0,x1,dx0,dx1,t0,t1,t)
    % Computes the Hermite Interpolant at t 
    % in a vector space.
    %   
    %   Input: Data x0, x1 in vector space format (ex. local coordinates)
    %          Derivative data dx0, dx1 in the same vector space
    %
    %   Output: Interpolated point at time t \in [t0,t1]
    %
    % 
    % a0 = @(t) 1 - 1/(t1-t0)^2*(t-t0)^2 + 2/(t1 - t0)^3*(t-t0)^2*(t-t1);
    % a1 = @(t) 1/(t1-t0)^2 * (t-t0)^2 - 2/(t1-t0)^3*(t-t0)^2*(t-t1);
    % b0 = @(t) (t-t0) - 1/(t1-t0)*(t-t0)^2 + 1/(t1-t0)*(t-t0)^2*(t-t1);
    % b1 = @(t) 1/(t1-t0)^2*(t-t0)^2*(t-t1);

    L00 = @(t) (1 - 2/(t0-t1)*(t-t0)) * ((t-t1)/(t0-t1))^2;
    L10 = @(t) (1 - 2/(t1-t0)*(t-t1)) * ((t-t0)/(t1-t0))^2;
    L01 = @(t) (t-t0)*((t-t1)/(t1-t0))^2;
    L11 = @(t) (t-t1)*((t-t0)/(t0-t1))^2;

    %gamma = @(t) a0(t)*x0 + a1(t)*x1 + b0(t)*dx0 + b1(t) * dx1;
    gamma = @(t) L00(t)*x0 + L10(t)*x1 + L01(t)*dx0 + L11(t) * dx1;
    y = gamma(t);
end