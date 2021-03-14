% This code assumes that CVX (http://cvxr.com/cvx/) is installed and available in the MATLAB path. 
% As default, the SDPT3 solver is used for CVX. 
% For our large scale setting, we observed that the Mosek solver (https://www.mosek.com/) for CVX usually performs better than the SDPT3 solver. 
% Please refer to http://cvxr.com/cvx/doc/mosek.html, which describes how to use Mosek solver in CVX. 
% Please refer to line 26 in proposed_train.m (proposed algorithm) and line 5 in compute_emd.m (EMD baseline) for changing the CVX solver option. 

clear
rng('default');
rseed = 1;
dimensions = 10;
M = 10;
N = M;

test_flag = true;
M_test = 200;
myRidge = 1e-8;
flag_squared = true;
sigma_rbf_vec = 0.1;
delta_factor = 100;

numdim = length(dimensions);
numseed = length(rseed);
length_m = length(M);
numsigma = length(sigma_rbf_vec);

for j = 1:numsigma
    sigma_rbf = sigma_rbf_vec(j);
    for d = 1:numdim
        d_src = dimensions(d);
        d_tgt = dimensions(d);
        for i = 1:length(M)
            m = M(i);
            n = N(i);
            for s = 1:numseed
                rng(rseed(s));
                fprintf('Synthetic data generation:\n');
                mu_src = zeros(d_src,1);
                mu_tgt = zeros(d_tgt,1);
                sigma_src = rand(d_src,d_src);
                sigma_src = sigma_src*sigma_src' + eye(d_src)*myRidge;
                sigma_src = (sigma_src+sigma_src')/2;
                sigma_src = sigma_src/trace(sigma_src);
                sigma_tgt = rand(d_tgt,d_tgt);
                sigma_tgt = sigma_tgt*sigma_tgt' + eye(d_tgt)*myRidge;
                sigma_tgt = (sigma_tgt+sigma_tgt')/2;
                sigma_tgt = sigma_tgt/trace(sigma_tgt);
                
                Xstr = mvnrnd(mu_src,sigma_src,m);
                Xttr = mvnrnd(mu_tgt,sigma_tgt,n);
                cost_Xstr_Xttr = euclidean_distances(Xstr,Xttr,flag_squared);
                kernel_Xstr = gaussianKernel(Xstr,Xstr,sigma_rbf);
                kernel_Xttr = gaussianKernel(Xttr,Xttr,sigma_rbf);
                if test_flag
                    Xste = mvnrnd(mu_src,sigma_src,M_test);
                    kernel_Xstr_Xste = gaussianKernel(Xstr,Xste,sigma_rbf);
                end
                %% Gaussian optimal - can be commented-out if not neeeded
                [opt_mapper] = gaussian_optimal(mu_src,mu_tgt,sigma_src,sigma_tgt);
                opt_pred = (mu_tgt + opt_mapper*( Xstr' - mu_src*ones(1,m) ))';
                %% EMD baseline - can be commented-out if not neeeded
                [alpha_emd, cost_emd] = emd_train(cost_Xstr_Xttr);
                emd_pred = barycenterSquaredEuclideanCost(alpha_emd,Xttr);
                emd_squared_error = sum((emd_pred-opt_pred).^2,2);
                emd_pred_error = mean(emd_squared_error);
                %% Proposed
                [cost_proposed, alpha_mat] = proposed_train(cost_Xstr_Xttr,kernel_Xstr,kernel_Xttr,delta_factor);                
                proposed_pred = barycenterSquaredEuclideanCost(alpha_mat,Xttr);
                proposed_squared_error = sum((proposed_pred-opt_pred).^2,2);
                proposed_pred_error = mean(proposed_squared_error);
                %% Out-of-sample (test) predictions
                if test_flag
                    opt_pred_test = mu_tgt + opt_mapper*( Xste' - mu_src*ones(1,M_test) );
                    opt_pred_test = opt_pred_test';% same orientation as Xste
                    pseudo_alpha_mat = kernel_Xstr_Xste'*(kernel_Xstr\alpha_mat);
                    proposed_pred_test = barycenterSquaredEuclideanCost(pseudo_alpha_mat,Xttr);
                    proposed_squared_error_test = sum((proposed_pred_test-opt_pred_test).^2,2);
                    proposed_pred_error_test = mean(proposed_squared_error_test);
                end
                fprintf('EMD mse: %g, Proposed mse: %g\n', emd_pred_error,proposed_pred_error);
                fprintf('Proposed out-of-sample mse: %g\n',proposed_pred_error_test);
            end
        end
    end
end
