function [mapper] = gaussian_optimal(mu_src,mu_tgt,sigma_src,sigma_tgt)
	%% compute optimal transport mapper for Gaussian-Gaussian case (Computational OT book)
	sigma_src_half = sqrtm(sigma_src);
	temp = sigma_src_half*sigma_tgt*sigma_src_half;
    temp = (temp + temp')/2;
	temp1 = sqrtm(temp);
	temp2 = inv(sigma_src_half);
	mapper = temp2*temp1*temp2;
end