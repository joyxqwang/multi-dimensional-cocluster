# multi_dimensional_cocluster

Code for the multi-dimensional co-clustering model.
This function optimizes the following problem:

    min_{F_i, G_i, S_i, F_ave} \sum_i ||X_i - F_i*S_i*G_i'||_1 +　\lambda * \sum_i||F_ave - F_i||_1
    
and

    min_{F_i, G_i, S_i, F_ave} \sum_i ||X_i - F_i*S_i*G_i'||_F^2 + \lambda * \sum_i||F_ave - F_i||_F^2

Format of input:

    n: number of samples
    d: number of SNPs
    c: number of QTs
    X: n*dim SNP data
    Y: n*c phenotype matrix
    options: a structure which optionally contain the following
        - N: numHidNodes
        - Act: choice of activation fucntion
        - computeAlpha: compute alpha if this value is 1
        - r: hyperparameter for the regularization term
                
Format of output:

    predFunc: A function handle to estimate the function for new points
    beta: weights of the hidden nodes, numHidNodes*c
    alpha: weights of each SNP in the prediction of phenotypes

Simply run the code in matlab as below:

    [predFunc, alpha] = AFNN(X, Y, options);
