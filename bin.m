

% Function to compute non-linear snapshots 
% as well as snapshots of the full system.





% use default tsteps/1000 for taking snapshots
N = 1024;

Tend = 8;
tsteps = 10e5;

% For time stepping 
dt = Tend/tsteps;

mys = 1:tsteps;
storeI = mys(mod(mys, tsteps/1000) == 0);


tList = 0:dt:Tend;


% Length 
L = 1;


% system parameters:
e = 0.015; 
b = 0.5;
gamma = 2;
c = 0.05;
h = L/(N + 2);


% Set up system matrix


% store only diagonal of E-matrix
E = ones(2*N,1);
% scale the upper part of E
E(1:N)= e*E(1:N);
E = (1/dt)*E; % For convenience when we do time-stepping


% Set up system matrix
K = toeplitz([-2 1 zeros(1, N - 2)]);
K(1, 1) = -1;
K(end, end) = -1;
K = -K;
A = [-e/h^2*K, -eye(N); b*eye(N), -gamma*eye(N)];
c = c*ones(2*N, 1);
A = sparse(A);
Yold = zeros(2*N, 1);
Y = [];
F = [];

nonlin = @(y) y.*(y-0.1).*(1-y);
bc = @(t) 50000*t^3*exp(-15*t);
for i=2:tsteps+1
    % Compute the nonlinear part 
    f = nonlin(Yold(1:N));
    
    g_plus_c_plus_f = c;
    
    g_plus_c_plus_f(1:N) = g_plus_c_plus_f(1:N) + f;
    
    g_plus_c_plus_f(1) = g_plus_c_plus_f(1) + e^2/h*bc(tList(i));
    
    rhs = A*Yold + g_plus_c_plus_f;
    Ynew = (1.0./E).*rhs + Yold;
    
    Yold = Ynew;
    if(find(storeI == i, 1))
        disp(tList(i));
        Y(:, end + 1) = Ynew;
        F(:,end + 1) = f; 

    end
end
% store the last snapshots
Y(:, end + 1) = Ynew;
F(:,end + 1) = f;    



plot(Y(1,:),Y(N+1,:))