function barycenters = barycenterSquaredEuclideanCost(alpha_weight,X)
	% alpha_weight: m x n
	% X: n x d
	n = size(alpha_weight,2);
	unnormalized_barycenters = alpha_weight*X;
	normalization = alpha_weight*ones(n,1);
	barycenters = bsxfun(@times,unnormalized_barycenters,1./normalization);
end

