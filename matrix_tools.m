% A structure containing some frequently used tools 
% in interpolation.

% Manifolds
%   Stiefel
%   Grassmann
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
    function D = dExpGrassmann(U,dir,dirtilde)
        % dir = Delta * t
        % dirtilde must be horizontal!

        [Q,S,V] = svd(dir,'econ');
        Scos = diag(cos(diag(S)));
        Ssin = diag(sin(diag(S)));

        %[Qd,Sd,Vd] = svd(dirtilde,'econ');
        
        [dQd,dSd,dVd] = dSVD2(Qd,Sd,Vd,dirtilde)

        Y = U * (V*Scos * V') + Q * (Ssin * V');

    end

MatTools.checkProjtan = @checkProjtangent;
    function e = checkProjtangent(U,v)
        P = U*U';
        w = (v*U'+U*v');
        e = norm(P*w+w*P-w);
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
MatTools.dSVD = @dSVD2;
function [Udot, Sdot, Vdot] = dSVDold(Y,r,dY)
    [U,S,V] = svd(Y,'econ');
    Ur = U(:,1:r); %U = Ur;
    S = diag(S);
    Sr = S(1:r); S = Sr;
    Vr = V(:,1:r); %V = Vr;
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
    for i = 1:r
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
    Gammar = 0.5 * (Gammar - Gammar');

    %Vdot = V * Gamma;
    
    Udot = (dY * Vr + Ur * (diag(Sr) * Gammar - Sdot)) / diag(Sr);

end
function [Udot, Sdot, Vdot] = dSVD(Y,r,dY)
%function [Udot, Sdot, Vdot] = dSVD(U,S,V,r,dY)
    [U,S,V] = svd(Y,'econ');
    
    Ur = U(:,1:r); %U = Ur;
    Sr = diag(S);
    Sr = Sr(1:r); S = Sr;
    Vr = V(:,1:r); %V = Vr;

    % Y = U * S * V', dY = Direction
    %
    % r is the number number determining the approximation Yr = UrSrVr'

    [n,~] = size(Ur);
    [m,~] = size(Vr);

    Sdot = zeros(r);
    for j = 1:r
        Sdot(j,j) = Ur(:,j)' * dY * Vr(:,j);
    end
    % Since Y most often will have low rank, one can make this 
    % computation more efficient. 
    Gamma = zeros(m,r);
    
    % The upper p x p block
    for i = 1:r
        w1 = Ur(:,i)'* dY ;
        w2 = dY * Vr(:,i);

        for j = 1:r
            if i ~= j
                if abs(S(i) - S(j)) < 10-10
                    error("The singular values for (i,j) = ("+num2str(Sr(i))+","+num2str(Sr(j))+") are too close")
                end
                Gamma(i,j) = Sr(i) * (Ur(:,i)'* dY * Vr(:,j)) + Sr(j) * (Ur(:,j)' * dY * Vr(:,i));

                Gamma(i,j) = Gamma(i,j) / ( (Sr(j) + Sr(i)) * (Sr(j) - Sr(i)) ) ;

            end

        end
    end

    Gammar = Gamma(1:r,1:r);
    Gammar = 0.5 * (Gammar - Gammar');

    %Vdot = V * Gamma;
    Vdot = 0;
    
    Udot = (dY * Vr + Ur * (diag(Sr) * Gammar - Sdot) ) / diag(Sr);

end

function [dU,dS,dV] = dSVD2(U,S,V,dY)
    % Y = U*S*V'
    [r,~] = size(S);
    [m,~] = size(V);
    S = diag(S);
    dS = zeros(r,1); 
    for i = 1:r
        dS(i) = U(:,i)' * dY * V(:,i);
    end
    dS = diag(dS);
    Gamma = zeros(m,r);
    for i = 1:m
            for j = 1:r
                if (i ~= j) && (i <= r)
                    if abs(S(i)-S(j))<10e-10
                        %error("The singular values are too close.")
                    end
                    Gamma(i,j) = S(i)*U(:,i)'*dY*V(:,j) + S(j)*U(:,j)'*dY*V(:,i);
                    Gamma(i,j) = Gamma(i,j) / ((S(j) + S(i))*(S(j)-S(i)));
                elseif i > r
                    Gamma(i,j) = U(:,j)'*dY*V(:,i) / S(j);
                end
            end
        
    end
    S = diag(S);
    dV = V * Gamma;
    dU = (dY * V(:,1:r) + U(:,1:r) * (S*Gamma(1:r,1:r) - dS)) / S;
    checkComp = 1;
    if checkComp
        % Check that the computation is correct
        dYsvd = dU * S * V(:,1:r)' + U(:,1:r) * dS * V(:,1:r)' + U(:,1:r) * S * dV';
        diff = norm(dY - dYsvd) / norm(dY);
        if diff < 10e-4
            fprintf("dSVD worked; difference %.2e \n",diff);
        else
            fprintf("dSVD did NOT work; difference %.2e \n",diff);
        end

    end
    
end

% MatTools.dSVD2 = @dSVD2;
% function [Udot, Sdot, Vdot] = dSVD2(Ur,Sr,Vr,r,dY)
%     % [U,S,V] = svd(Y,'econ');
%     % Ur = U(:,1:r);
%     % S = diag(S);
%     % Sr = S(1:r);
%     % Vr = V(:,1:r);
%     % Yr = Ur * Sr * Vr', dY = Direction
%     %
%     % r is the number number determining the approximation Yr = UrSrVr'
% 
%     [n,~] = size(Ur);
%     [m,r] = size(Vr);
% 
%     Sdot = zeros(r);
%     for j = 1:r
%         Sdot(j,j) = Ur(:,j)' * dY * Vr(:,j);
%     end
%     % Since Y most often will have low rank, one can make this 
%     % computation more efficient. 
%     Gamma = zeros(m,r);
%     for i = 1:m
%         w1 = U(:,i)'* dY ;
%         w2 = dY * V(:,i);
% 
%         for j = 1:r
%             if i ~= j
%                 if abs(S(i) - S(j)) < 10-10
%                     error("The singular values for (i,j) = ("+num2str(Sr(i))+","+num2str(Sr(j))+") are too close")
%                 end
%                 Gamma(i,j) = S(i) * w1 * V(:,j) + S(j) * U(:,j)' * w2;
% 
%                 Gamma(i,j) = Gamma(i,j) / ( (S(j) + S(i)) * (S(j) - S(i)) ) ;
%             end
% 
%         end
%     end
%     Gammar = Gamma(1:r,1:r);
%     Gammar = 0.5 * (Gammar - Gammar');
% 
%     Vdot = V * Gamma;
% 
%     Udot = (dY * Vr + Ur * (diag(Sr) * Gammar - Sdot)) / diag(Sr);
% 
% end
% 
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



%% Aux. matrix tools
MatTools.house_qr = @house_qr;
function [U, R, Q] = house_qr(X)
% Stewart, Mat Decomp 1, p. 259
%
%
[n,k] = size(X);

% store Householder vectors
U = zeros(n,k);
R = zeros(k,k);
Q = eye(n);
for j=1:k
    % cancel column j
    [U(j:n,j), R(j,j)] = housegen(X(j:n,j));
    vT                 = U(j:n,j)'*X(j:n, j+1:k);
    X(j:n,j+1:k)       = X(j:n,j+1:k) - U(j:n,j)*vT;
    R(j,j+1:k)         = X(j,j+1:k);
    %
    % tmp
    %Q = Q*(eye(n) - U(:,j)*U(:,j)');
end
end

MatTools.housegen = @housegen;
function [u,nu] = housegen(x)
% From Stewart: Mat Decomp1, p. 257
%This algorithm takes a vector x and produces a vector u 
% that generates a Householder transformation H = I - uu^T 
% such that Hx = = +/-||x'|| ei. The quantity +/- ||x|| is
%returned in v
%
%
u = x;
nu= norm(u, 2);
if nu == 0
    u(l) = sqrt(0.5);
    return;
end

u = x/nu;
if (u(1)>0)
    u(1) = u(1) + 1;
    nu = -nu;
else
    u(1) = u(1)-1;
end

u = u/(sqrt(abs(u(1))));

end

MatTools.house_block = @house_block;
function [W,Y] = house_block(U)
% low rank representation of Q-factor of QR  
% from Householder approach
%
% implementation follows Golub/Van Loan/4th edition,
% Algorithm 5.1.2, p. 239
%

[n,k] = size(U);

W = zeros(n,k);
Y = zeros(n,k);

W(:,1) = U(:,1);
Y(:,1) =  U(:,1);
for j=2:k
    % compute z = -Q*U(:,j) = -(I+WY')*U(:,j)
    z      = (U(:,j)-W(:,1:j-1)*(Y(:,1:j-1)'*U(:,j)));
    W(:,j) = z;
    Y(:,j) = U(:,j);
end

end

MatTools.applyWYT = @appplyWYT;
function QU = appplyWYT(U,W,Y)
    QU = U - Y*(W'*U);
end

MatTools.applyWY = @appplyWY;
function QU = appplyWY(U,W,Y)
    QU = U - W*(Y'*U);
end

MatTools.findOptQ = @optimalQ;
function [W,Y] = optimalQ(U1,U2,W1,Y1,W2,Y2)
    % Find the optimal Q = I + YW'.
    Q1U1 = appplyWYT(U1,W1,Y1);
    Q1U2 = appplyWYT(U2,W1,Y1);
    Q2U1 = appplyWYT(U1,W2,Y2);
    Q2U2 = appplyWYT(U2,W2,Y2);

    [~,p] = size(U1);
    for i = 1:4
         
    end
end

end