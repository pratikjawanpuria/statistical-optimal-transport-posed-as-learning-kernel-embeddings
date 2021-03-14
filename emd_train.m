function [coupling, cost_emd] = emd_train(cost)
    M = size(cost,1);
    N = size(cost,2);
    a = (1/M)*ones(M,1);
    b = (1/N)*ones(N,1);
    coupling = compute_emd(a, b, cost);
    cost_emd = coupling(:)'*cost(:);
end

