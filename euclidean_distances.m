function cost = euclidean_distances(X,Z,flag_squared)
    M = size(X,1);
    N = size(Z,1);
    XX = sum(X.*X,2)*ones(1,N);
    ZZt = (sum(Z.*Z,2)*ones(1,M))';
    distance = -2*(X*Z') + XX + ZZt;
    if flag_squared
        cost = distance;
    else
        cost =sqrt(distance);
    end
end

