function kernel_matrix = get_kernel_distance(kernel_X, kernel_Z, kernel_XZ, flag_squared)
    a = diag(kernel_X);
    b = diag(kernel_Z);
    M = size(a,1);
    N = size(b,1);
    kernel_matrix = a*ones(1,N) + (b*ones(1,M))' - 2*kernel_XZ;
    if flag_squared
    else
        kernel_matrix = sqrt(kernel_matrix);
    end
end