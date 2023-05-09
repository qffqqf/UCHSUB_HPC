function [M, K] = clean_matrices(M, K, threshold, dofs_red)

% M(dofs_red.I,dofs_red.I) = take_diag(M(dofs_red.I,dofs_red.I));
% K(dofs_red.I,dofs_red.I) = take_diag(K(dofs_red.I,dofs_red.I));
% K(dofs_red.I,dofs_red.A) = 0;
% K(dofs_red.A,dofs_red.I) = 0;
M(dofs_red.I,dofs_red.I) = clean_zeros(M(dofs_red.I,dofs_red.I),1e10);
K(dofs_red.I,dofs_red.I) = clean_zeros(K(dofs_red.I,dofs_red.I),1e10);

M = clean_zeros(M,threshold);
K = clean_zeros(K,threshold);