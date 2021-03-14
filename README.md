This package contains a MATLAB implementation of the kernel embedding based optimal transport algorithm proposed in [1]. 

[1] J. Saketha Nath and Pratik Jawanpuria. Statistical Optimal Transport posed as Learning Kernel Embedding. In Conference on Neural Information Processing Systems (NeurIPS), 2020.

Please cite the above paper if you are using this code.

(c) 2020-21 J. Saketha Nath <saketha@cse.iith.ac.in> and Pratik Jawanpuria <pratik.iitb@gmail.com>


Usage of code: 
--------------
The file test.m creates a synthetic dataset (as described in [1]) and computes the results with the proposed algorithm and the earth mover distance baseline. Our implementation assumes that CVX (http://cvxr.com/cvx/) is installed and available in the MATLAB path. 
