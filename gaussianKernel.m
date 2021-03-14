function K = gaussianKernel(X,Z,sigma)
    flag_squared = true;
    K = exp(-euclidean_distances(X,Z,flag_squared)/(2*sigma*sigma));
end