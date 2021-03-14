function coupling = compute_emd(a,b,cost)
    M = size(a,1);
    N = size(b,1);
    cvx_begin quiet 
        cvx_solver sdpt3 %mosek
        variable P(M,N)
        minimize sum(sum(P.*cost))
        subject to
            P*ones(N,1)==a;
            P'*ones(M,1)==b;
            P>=0;
    cvx_end
    coupling = P;
end