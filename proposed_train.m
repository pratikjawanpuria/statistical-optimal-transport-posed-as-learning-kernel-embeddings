function [cost_proposed, alpha_var_opt] = proposed_train(cost_Xstr_Xttr,kernel_Xstr,kernel_Xttr,delta_factor)
    m = size(cost_Xstr_Xttr,1);
    n = size(cost_Xstr_Xttr,2);
    ones_m = ones(m,1);
    ones_n = ones(n,1);
    
    kernel_Xstr_sqrt = sqrtm(kernel_Xstr);
    kernel_Xstr_sqrt = (kernel_Xstr_sqrt+kernel_Xstr_sqrt')/2;
    kernel_Xttr_sqrt = sqrtm(kernel_Xttr);
    kernel_Xttr_sqrt = (kernel_Xttr_sqrt+kernel_Xttr_sqrt')/2;

    Hkernel_Xstr_sqrt = sqrtm(kernel_Xstr.*kernel_Xstr);
    Hkernel_Xstr_sqrt = (Hkernel_Xstr_sqrt+Hkernel_Xstr_sqrt')/2;
    Hkernel_Xttr_sqrt = sqrtm(kernel_Xttr.*kernel_Xttr);
    Hkernel_Xttr_sqrt = (Hkernel_Xttr_sqrt+Hkernel_Xttr_sqrt')/2;

    myeps_m = (sum(sum(kernel_Xstr))/m^2)/delta_factor;
    myeps_n = (sum(sum(kernel_Xttr))/n^2)/delta_factor;
    Hmyeps_m = (sum(sum(kernel_Xstr.*kernel_Xstr))/m^2)/delta_factor;
    Hmyeps_n = (sum(sum(kernel_Xttr.*kernel_Xttr))/n^2)/delta_factor;
    
%     G1inv = inv(kernel_Xstr);
%     G2inv = inv(kernel_Xttr);
    cvx_begin
        cvx_precision low
        cvx_solver sdpt3 %mosek
        variable alpha_var(m,n);
        
        expression b1(m,1);
        expression b2(n,1);
        expression a1(m,n);
        expression a2(m,n);

        b1 = (alpha_var*ones(n,1) - (1/m)*ones_m);
        b2 = (alpha_var'*ones(m,1) - (1/n)*ones_n);
        
%         a1 = G1inv*alpha_var;
%         a2 = alpha_var*G2inv;
        a1 = kernel_Xstr\alpha_var;
        a2 = alpha_var/kernel_Xttr;

        minimize trace(cost_Xstr_Xttr'*alpha_var)

        subject to
            norm( kernel_Xstr_sqrt * b1 ) <= sqrt(myeps_m);
            norm( kernel_Xttr_sqrt * b2 ) <= sqrt(myeps_n);
            %quad_form(b1,kernel_Xstr) <= myeps_m;
            %quad_form(b2,kernel_Xttr) <= myeps_n;
            norm( Hkernel_Xstr_sqrt * b1 ) <= sqrt(Hmyeps_m);
            norm( Hkernel_Xttr_sqrt * b2 ) <= sqrt(Hmyeps_n);

            a1 >= zeros(m,n);
            a2 >= zeros(m,n);

            sum(sum(alpha_var)) == 1;
            alpha_var >= zeros(m,n);
    cvx_end

    cost_proposed = cvx_optval;
    alpha_var_opt = alpha_var;
end
    
    