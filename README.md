# multi_dimensional_cocluster

Code for the multi-dimensional co-clustering model.
These functions optimize the following problem:

1) ccl_L1.m:

    min_{F_i, G_i, S_i, F_ave} \sum_i ||X_i - F_i*S_i*G_i'||_1 +　\lambda * \sum_i||F_ave - F_i||_1
    
and

2) ccl_L2.m:

    min_{F_i, G_i, S_i, F_ave} \sum_i ||X_i - F_i*S_i*G_i'||_F^2 + \lambda * \sum_i||F_ave - F_i||_F^2

Format of input:

    n: number of samples
    d: number of features
    X: cell array, each array denotes a d*n_i matrix for the i-th data type
    NClusters: number of cross-tissue cancer clusters
    NClustersG: number of feature clusters
    iniF: initialization of cross-tissue cluster indicator
    lambda: initial value of the hyper-parameter, default value is 10^5

Format of output:

    outF: cross-tissue cancer cluster indicator from each data type
    outG: feature cluster indicator of each data type
    outF_ave: integrated cross-tissue cancer cluster indicator
    obj: objective function value
    
To use the function, simply run the code in matlab as below:

    [outF, outG, S, outF_ave, obj] = ccl_L1(X, NClusters, NClustersG, iniF)

or

    [outF, outG, S, outF_ave, obj] = ccl_L2(X, NClusters, NClustersG, iniF)
