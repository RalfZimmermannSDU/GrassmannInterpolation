% A structure containing some frequently used tools 
% in interpolation.

% Manifolds
%
%   Stiefel S
%   Grassmann G
%   Symplectic Stiefel Sp
%
%
function MatTools = matrix_tools()

%% Stiefel

% Generate Stiefel matrix
MatTools.RandS = @RandStiefel;
    function Q = RandStiefel(n,p)
        [Q,~] = qr(rand(n,p),'econ');
    end


% Generate Stiefel tangent at a point Q
MatTools.vectorS = @vectorStiefel;
    function v = vectorStiefel(U)
        [n,p] = size(U);
        A = rand(p,p);
        A = 0.5 * (A-A');
        T = rand(n,p);
       
        v = U*A + (eye(n) - U*U')*T; % Efficiency does not matter here
    end

% Stiefel geodesics
MatTools.ExpS = @ExpStiefel;
    function Y = ExpStiefel(U,Dir,t)
        if nargin < 3 || isempty(t)
            t = 1;
        end
        [~,p] = size(U);
        
        A = U'*Dir;
        UcomplT = Dir - U*A;
        [Q,R] = qr(UcomplT,'econ');
    
        Y = [U Q] * expm(t * [A -R'; R zeros(p)]) * [eye(p); zeros(p,p)];
    end

% Stiefel logarithm
MatTools.LogS = @LogStiefel;
    function Log = LogStiefel(U,tildeU,tau)
        % Stiefel logarithm with anchor U and input V
        if nargin < 3
            tau = 1e-13;
        end
        maxiter = 30;
    
        [~,p] = size(U);
        M = U'*tildeU;
        [Q,N] = qr(tildeU - U * M,'econ');
        
        [Vk,~] = qr([M; N]);
        if norm(Vk(:,1:p) - [M;N],'fro') > 1e-11
            Vk = -Vk; % Fix sign change caused by QR
            if norm(Vk(:,1:p) - [M;N],'fro') > 1e-11
                error("Error in initializing Stiefel Log ")
            end
        end
        % Fix via Procustes of Y0 
        for k = 1:maxiter
            L = logm(Vk);
            if norm(L(p+1:2*p,p+1:2*p)) < tau
                Log = L;
                break;
            end
        end
    end

% Check if a matrix is in the Stiefel manifold
MatTools.CheckS = @CheckStiefel;
    function CheckStiefel(U)
        [~,p] = size(U);
    
        w = norm(U'*U - eye(p),'fro');
        fprintf("||U'U - I|| = %2.3e \n",w);
    end


%% Grassmann 

% Generate Grassmann matrix
MatTools.RandG = @RandGrassmann;
    function U = RandGrassmann(n,p)
        [U,~] = qr(rand(n,p),'econ');
    end

% Generate Grassmann tangent at U
MatTools.vectorG = @vectorGrassmann;
    function v = vectorGrassmann(U)
        % We want to work with n x p - matrices, so we use the orthogonal
        % projection via lifts to the Stiefel manifold
        Z = rand(size(U));
        v = Z - U * (U'*Z);
    end

% Grassmann geodesics
MatTools.ExpG = @ExpGrassmann;
    function Y = ExpGrassmann(U,dir,t)
        % This is equation 3.10 in Bendokat and Zimmermann 2024
        if nargin < 3 || isempty(t)
            t = 1;
        end
        Dir = t * dir;
        [Q,S,V] = svd(Dir,'econ');
        Scos = diag(cos(diag(S)));
        Ssin = diag(sin(diag(S)));

        Y = U * (V*Scos * V') + Q * (Ssin * V');

        %[Y,~] = qr(Y,'econ');
    end

MatTools.dExpG = @dExpGrassmann;
    function D = dExpGrassmann()

    end

% Grassmann log
MatTools.LogG = @LogGrassmann;
    function Log = LogGrassmann(U,Y)
        % This is Algorithm 1 in Bendokat and Zimmermann 2024
        [Q,~,R] = svd(Y'*U);
        Ystar = Y * (Q*R');
        
        IUUTY = Ystar - U * (U'*Ystar);

        [Q,S,R] = svd(IUUTY,'econ');

        Sigma = diag(asin(diag(S)));

        Log = Q*Sigma*R';
    end

% Local cooridnates and parametrization
MatTools.LocalCoordG = @LocalCoordinatesGrassmann;
    function [BAinv,Acond] = LocalCoordinatesGrassmann(P,n,p)
        %
        % Compute the local coordinates of a point on Gr(n,p)
        %
        % P = (A B')
        %     (B C )
        %
        % P = U*U';
        %
        % OR 
        % Inserting a Stiefel representative U st. P = U*U' 
        % 
        
        A = P(1:p,1:p);
        B = P(p+1:n,1:p);

        BAinv = B / A;
        Acond = cond(A);
    end

MatTools.dLocalCoordG = @dLocalCoordinatesGrassmann;
    function [dBAinv] = dLocalCoordinatesGrassmann(P,Delta,n,p)
        %
        % Compute the derivaitve of the local coordinates 
        % of a curve on Gr(n,p)
        %
        % P = (A B')
        %     (B C )
        %
        % P = U*U' is the base point
        % Delta is the direction

        PA = P(1:p,1:p);
        PB = P(p+1:n,1:p);

        DA = Delta(1:p,1:p);
        DB = Delta(p+1:n,1:p);
        
        dBAinv = DB / PA - PB * (PA \ (DA / PA));
    end

MatTools.ParamG = @ParametrizationGrassmann;
    function P = ParametrizationGrassmann(BAinv)
        % Compute the projector P from the local coordinates BAinv
        [~,p] = size(BAinv);

        P = ([eye(p); BAinv] / (eye(p) + (BAinv)'*BAinv) ) * [eye(p) BAinv'];
    end

MatTools.dParamG = @dParametrizationGrassmann;
    function dP = dParametrizationGrassmann(B,Delta)
        % Compute the derivative of the parametrization
        % Follows the computaton provided in the shared document
     
        [nmp,p] = size(B);
        n = nmp + p;

        IBBT = eye(p) + B'*B;
        P = ParametrizationGrassmann(B);
        DinvIBBT = Delta / IBBT;
        DinvIBBTB = DinvIBBT*B';

        M1 = zeros(n,n);
        M2 = zeros(n,n);

        M1(p+1:end,1:p) = DinvIBBT;
        M1(p+1:end,p+1:end) = DinvIBBTB;
        
        M2(1:p,p+1:end) = Delta';
        M2(p+1:end,1:p) = Delta;
        
        M2 = P * M2 * P;

        dP = M1 - M2 + M1';
    end

%MatTools.dSVD = @dSVD;
% Differentiate the truncated singular value decomposition
MatTools.dSVD = @dSVD;
function [Udot, Sdot, Vdot] = dSVD(Y,r,dY)
    [U,S,V] = svd(Y);
    Ur = U(:,1:r);
    S = diag(S);
    Sr = S(1:r);
    Vr = V(:,1:r);
    % Y = U * S * V', dY = Direction
    %
    % r is the number number determining the approximation Yr = UrSrVr'

    [n,~] = size(Ur);
    [m,r] = size(Vr);
    
    Sdot = zeros(r);
    for j = 1:r
        Sdot(j,j) = Ur(:,j)' * dY * Vr(:,j);
    end
    
    % Since Y most often will have low rank, one can make this 
    % computation more efficient. 
    Gamma = zeros(m,r);
    for i = 1:m
        w1 = U(:,i)'* dY ;
        w2 = dY * V(:,i);
        
        for j = 1:r
            if i ~= j
                if abs(S(i) - S(j)) < 10-10
                    error("The singular values for (i,j) = ("+num2str(Sr(i))+","+num2str(Sr(j))+") are too close")
                end
                Gamma(i,j) = S(i) * w1 * V(:,j) + S(j) * U(:,j)' * w2;

                Gamma(i,j) = Gamma(i,j) / ( (S(j) + S(i)) * (S(j) - S(i)) ) ;
            end

        end
    end
    Gammar = Gamma(1:r,1:r);

    Vdot = V * Gamma;
    
    Udot = (dY * Vr + Ur * (diag(Sr) * Gammar - Sdot)) / diag(Sr);

end

MatTools.dQR = @dQR;
function [dQ,dR] = dQR(Y,dY)
    [n,p] = size(Y);
    
    [Q,R] = qr(Y,'econ');
    
    PL = tril(ones(p,p),-1);
    L = PL.*(Q'*dY / R);
    X = L - L';
    
    dR = Q'*dY - X*R;
    dQ = (eye(n) - Q*Q') * dY / R + Q*X;

end

MatTools.LC_distbound = @lcdistbound;
function cb = lcdistbound(U,V)
    % U and V are Stiefel matrices
    [~,p] = size(U);
    U1 = U(1:p,1:p);
    V1 = V(1:p,1:p);

    cb = sqrt(p*norm(inv(U1),'fro')^2-cond(U1,'fro')^2)+sqrt(p*norm(inv(V1),'fro')^2-cond(V1,'fro')^2);
    
end


MatTools.dphi_cond_bound = @conditionnumber_phi_bound;
function cb = conditionnumber_phi_bound(B)
    S = svd(B,'econ');
    mx = max(S.^2 ./ ((1+S.^2).^2));
    
    cb = sqrt(2)*sqrt( 1/(1+S(end)^2)^2 + mx )+1;
end





end