% min_{F_i, G_i, S_i, F_ave} \sum_i ||X_i - F_i*S_i*G_i'||_F^2 +
% \lambda * \sum_i||F_ave - F_i||_F^2
function [outF, outG, S, outF_ave, obj, lambda] = ccl_L2(X, NClusters, NClustersG, iniF, lambda)
% X: contains nType different type of d*n data
% F_i: row cluster indicator matrix, d*NClusters
% G_i: column cluster indicator matrix, n*NClusters
% F_ave: average of the row cluster indicators
% Time complexity: O(NdcT), where N = \sum_i n_i

% path_data = 'F:\JoyWang\UTA\PhD\TCGA Landscape\TCGA_landscape_data\';
% dataList = {'BRCA', 'COAD', 'GBM', 'KIRC', 'KIRP', 'LGG', 'LUAD', 'LUSC', 'OV', 'READ', 'UCEC'};

% Initialization
NITER = 100;
nType = length(X);
d = size(X{1}, 1);
F_ave = zeros(d, NClusters);
if nargin < 4
    lambda = 1e5;
end

for ii = 1 : nType
    % Initialize F and G
    XX = X{ii};
    n = size(XX, 2);
    nRepeat = 100;
	ncG = NClustersG(ii);
    iniG = zeros(ncG, nRepeat);
	%iniG = zeros(n, nRepeat);
	for rr = 1 : nRepeat
		%iniG(:, rr) = randsrc(n, 1, 1:NClusters);
        tmp = randperm(n);
        iniG(:, rr) = tmp(1:ncG);
	end;
%     [bestIndF] = tuneKmeans(XX, iniF);
%     [bestIndG] = tuneKmeans(XX', iniG);
    [bestIndF] = tuneKmeans_new(XX, iniF);
    [bestIndG] = tuneKmeans_new(XX', iniG);
    F_init = [];
    G_init = [];
    for jj = 1 : NClusters
        F_init(bestIndF==jj,jj) = 1;
    end
    for jj = 1 : ncG
        G_init(bestIndG==jj,jj) = 1;
    end
    F{ii} = F_init;
    G{ii} = G_init;
    invF = diag(1./(sum(F{ii},1)+eps));
    invG = diag(1./(sum(G{ii},1)+eps));
        
    % Initialize S
    for kk = 1 : NClusters
        tmpF(kk, :) = sum(XX(F{ii}(:,kk) == 1, :), 1);
    end
    for kk = 1 : ncG
        tmpG(:, kk) = sum(tmpF(:, G{ii}(:,kk) == 1), 2);
    end
    S{ii} = invF * tmpG * invG;
    clear tmpF tmpG;
    
    tmpLoss = norm(XX - F{ii}*S{ii}*G{ii}','fro');
    tmpReg = norm(F{ii}-F_ave, 'fro');
    tmpObj(ii) = tmpLoss^2 + lambda*tmpReg^2;
end
% Initialize F_ave
for ii = 1 : nType
    F_ave = F_ave + F{ii};
end
F_ave = F_ave / nType;

% Compute Obj
obj(1) = sum(tmpObj);

for it = 2 : NITER
    
    for ii = 1 : nType
        % Initialize F and G
        XX = X{ii};
        n = size(XX, 2);
        ncG = NClustersG(ii);
        for jj = 1 : d
            idx = solveF(XX(jj,:), S{ii}*G{ii}', F_ave(jj,:), lambda);
            tmp = zeros(1, NClusters);
            tmp(idx) = 1;
            F{ii}(jj,:) = tmp;
        end
        for jj = 1 : n
            idx = solveG(XX(:,jj), F{ii}*S{ii});
            tmp = zeros(1, ncG);
            tmp(idx) = 1;
            G{ii}(jj,:) = tmp;
        end
        invF = diag(1./(sum(F{ii},1)+eps));
        invG = diag(1./(sum(G{ii},1)+eps));
        
        % Initialize S
        for kk = 1 : NClusters
            tmpF(kk, :) = sum(XX(F{ii}(:,kk) == 1, :), 1);
        end
        for kk = 1 : ncG
            tmpG(:, kk) = sum(tmpF(:, G{ii}(:,kk) == 1), 2);
        end
        S{ii} = invF * tmpG * invG;
        clear tmpF tmpG;
        tmpLoss = norm(XX - F{ii}*S{ii}*G{ii}','fro');
        tmpReg = norm(F{ii}-F_ave, 'fro');
        tmpObj(ii) = tmpLoss^2 + lambda*tmpReg^2;
    end
    
    % Update F_ave
    F_ave = zeros(d, NClusters);
    for ii = 1 : nType
        F_ave = F_ave + F{ii};
    end
    F_ave = F_ave / nType;
    if isempty(find(F_ave>0 & F_ave<1)) == 0
        lambda = lambda * 2;
    end
    
    % Update Obj
    obj(it) = sum(tmpObj);
    if mod(it, 10) == 0
        fprintf('%dth iteration, obj = %f \n', it, obj(it));
    end
    
    if abs(obj(it) - obj(it-1)) < 10^-6*obj(it-1)
        break;
    end
    
end

%% Output
for ii = 1 : nType
    ncG = NClustersG(ii);
    for jj = 1 : NClusters
        tmpF(F{ii}(:,jj)==1) = jj;
    end
    for jj = 1 : ncG
        tmpG(G{ii}(:,jj)==1) = jj;
    end
    outF{ii} = tmpF';
    outG{ii} = tmpG';
    clear tmpF;
    clear tmpG;
end
for jj = 1 : NClusters
    outF_ave(F_ave(:,jj)==1,1) = jj;
end

end