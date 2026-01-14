%% Interpolate the FN system
%
% Interpolate the POD basis of the u-component of the FN system
%   Data: Contains U 
%   dData: Contains the dU's, which are all horizontal.
%   Data_ref: Contains reference U data, matching the U for each m_h entry
%               eg. Data_ref{m_k*k + 1).U = Data{k}.U
function FN_interpolate_update(Data, dData, Data_ref, points, m_h, p, maxsteps)

M = matrix_tools();

E_lag_loc = [];
E_herm_loc = [];
E_lag_norm = [];
E_herm_norm = [];

% Interpolate piece-wise on points(j) -- points(j+1)

l = length(points);


% Normal coordinates
for i = 1:(l-1)
    ts = linspace(points(i),points(i+1),m_h+1);
    Ds = Data(i:(i+1));
    dDs = dData(i:(i+1));
    Dref = Data_ref((i-1)*m_h + 1:m_h*i+1);
    for j = 1:m_h+1
        Ptrue = Dref{j}(:,1:p)*Dref{j}(:,1:p)';
        Uint = Interpolate_Gr([points(i) points(i+1)],Ds,ts(j),'normal_lag');
        e = norm(Uint*Uint' - Ptrue,'fro') / norm(Ptrue,'fro');
        E_lag_norm = [E_lag_norm e];
    end
    if i ~= l-1
        E_lag_norm = E_lag_norm(1:end-1);
    end
end

for i = 1:(l-1)
    ts = linspace(points(i),points(i+1),m_h+1);
    Ds = Data(i:(i+1));
    dDs = dData(i:(i+1));
    Dref = Data_ref((i-1)*m_h + 1:m_h*i+1);
    for j = 1:m_h+1
        Ptrue = Dref{j}(:,1:p)*Dref{j}(:,1:p)';
        Uint = Interpolate_Gr([points(i) points(i+1)],Ds,ts(j),'normal_herm',dDs);
        e = norm(Uint*Uint' - Ptrue,'fro') / norm(Ptrue,'fro');
        E_herm_norm = [E_herm_norm e];
    end
    if i ~= l-1
        E_herm_norm = E_herm_norm(1:end-1);
    end
end


% Maxvol coordinates
% First we compute the maximum-volume 
% One has to selet the permutation matrix resulting in the smallest maximal
% conditon number, for each pair of consecutive data points. 
% For this example, it is seen that it makes no large difference if one
% chooses the P obtained at entry s = 1 or s = 2. Choose s = 1.
conds = zeros(l,4);
P = cell(1,l-1);
W = cell(1,l-1);
Y = cell(1,l-1);
for i = 1:(l-1)
    ss = i:(i+1);
    Ps = cell(1,2);
    for s = 1:2
        [~,Ps{s}] = maxvol(Data{ss(s)},maxsteps);
        %conds(i,s) = cond(Data{ss(s)}(1:p,1:p),'fro');
        PU = Ps{s}*Data{ss(1)};
        conds(i,(s-1)*2+1) = cond(PU(1:p,1:p) ,'fro');
        PU = Ps{s}*Data{ss(2)};
        conds(i,(s-1)*2+2) = cond(PU(1:p,1:p),'fro');
    end
    %P{i} = Ps{1};
    %P{i} = Ps{2};
    [UQ,R,Q] = M.house_qr(Data{ss(2)});
    [W{i},Y{i}] = M.house_block(UQ);

end


for i = 1:(l-1)
    ts = linspace(points(i),points(i+1),m_h+1);
    Ds = Data(i:(i+1));
    dDs = dData(i:(i+1));

    % Apply P
    for s = 1:2
        %Ds{s} = P{i}*Ds{s};
        %dDs{s} = P{i}*dDs{s};
        Ds{s} = M.applyWYT(Ds{s},W{i},Y{i});
        dDs{s} = M.applyWYT(dDs{s},W{i},Y{i});
    end
    Dref = Data_ref((i-1)*m_h + 1:m_h*i+1);
    for j = 1:m_h+1
        Ptrue = Dref{j}(:,1:p)*Dref{j}(:,1:p)';
        Uint = Interpolate_Gr([points(i) points(i+1)],Ds,ts(j),'local_lag');
        % Uint = P{i}'*Uint;
        Uint = M.applyWY(Uint,W{i},Y{i});
        e = norm(Uint*Uint' - Ptrue,'fro') / norm(Ptrue,'fro');
        E_lag_loc = [E_lag_loc e];
    end
    if i ~= l-1
        E_lag_loc = E_lag_loc(1:end-1);
    end
end

for i = 1:(l-1)
    ts = linspace(points(i),points(i+1),m_h+1);
    Ds = Data(i:(i+1));
    dDs = dData(i:(i+1));
    norm(Ds{s}'*dDs{s});
    % Apply P
    for s = 1:2
        % Ds{s} = P{i}*Ds{s};
        % dDs{s} = P{i}*dDs{s};
        Ds{s} = M.applyWYT(Ds{s},W{i},Y{i});
        dDs{s} = M.applyWYT(dDs{s},W{i},Y{i});
    end
    Dref = Data_ref((i-1)*m_h + 1:m_h*i+1);
    for j = 1:m_h+1
        Ptrue = Dref{j}(:,1:p)*Dref{j}(:,1:p)';
        Uint = Interpolate_Gr([points(i) points(i+1)],Ds,ts(j),'local_herm',dDs);
        %Uint = P{i}'*Uint;
        Uint = M.applyWY(Uint,W{i},Y{i});
        e = norm(Uint*Uint' - Ptrue,'fro') / norm(Ptrue,'fro');
        E_herm_loc = [E_herm_loc e];
    end
    if i ~= l-1
        E_herm_loc = E_herm_loc(1:end-1);
    end
end


figure
plot(E_lag_loc)
hold on
plot(E_lag_norm)
hold on
plot(E_herm_loc)
plot(E_herm_norm)
legend('Lag MV','Lag RN','Herm MV','Lag RN')

end