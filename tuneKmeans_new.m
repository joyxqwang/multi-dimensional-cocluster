% M: n*d data matrix for clustering
% for function "litekmeans" please refer to http://www.cad.zju.edu.cn/home/dengcai/Data/Clustering.html
function [finalInd Ind kmobj min] = tuneKmeans_new(M, Ini)

min = Inf;
nIni = size(Ini, 2);
NClusters = size(Ini, 1);
kmobj = zeros(nIni);

for ii = 1 : nIni
    [Ind(:,ii), center, bCon, sumD, D] = litekmeans(M, NClusters, 'Start', M(Ini(:,ii), :));
    obj = sum(sumD);
    kmobj(ii) = obj;
    if obj < min
        min = obj;
        finalInd = Ind(:, ii);
    end;
end;

end

