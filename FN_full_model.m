% Simulate the FitzHigh-Nagumo model
function Y = FN_full_model(Ia)

if nargin < 1
    Ia = 0.03;
end

T = 8; % Final time
L = 1; % x \in [0,L]

nt = 10e5;
nx = 1024;

dx = L /(nx + 2); % spatial stepsize
dt = T / nt; % time stepsize

% Constants
epsilon = 0.015;
b = 0.3;
gamma = 0.9;
%Ia = 0.05;


% Set up system matrix 
M = 2*diag(ones(1,nx))-diag(ones(1,nx-1),1)-diag(ones(1,nx-1),-1);
M(1,1)=1; M(end,end) = 1;
M = -epsilon^2 * 1/(dx^2)*M; 

A = [M, -eye(nx); b*eye(nx), -gamma*eye(nx)];
A = sparse(A);
Ia = Ia*ones(2*nx, 1);
beta = @(t) 50000 * t^3 * exp(-15*t);

% Initial value
y = zeros(2*nx,1);

Y = [y]; % Snapshots
F = [y]; % Snapshots nonlinear (if need be)

% take snapshots for at each 1000 steps, i.e. at k*1000 == 0 mod(1000)
StoreIndices = 1:nt;
StoreIndices = StoreIndices(mod(StoreIndices,nt/1000) == 0);


for i = 1:nt
    t = i*dt;

    f = y(1:nx).*(y(1:nx) - 0.1).*(1 - y(1:nx));

    g = Ia; 
    g(1:nx) = g(1:nx) + f;
    g(nx+1:end) = zeros(nx,1);

    g(1) = g(1) + epsilon^2 / dx * beta(t); % Boundary cond. in pos. 1

    rhs = A*y + g;
    rhs(1:nx) = 1/epsilon * rhs(1:nx);
    y = y + rhs * dt;
    if find(StoreIndices == i, 1)
        %disp(t);
        Y(:,end+1) = y;
    end
end




% figure
% for j = 1:3:nx
%     plot(Y(j,:),Y(j + nx,:),'r');
%     hold on
% end

end