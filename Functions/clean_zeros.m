function M = clean_zeros(M, thres)

scalar = max(max(M));
% id_delete = find(abs(M) < scalar/thres);
M(abs(M) < scalar/thres) = 0;
M = sparse(M);
