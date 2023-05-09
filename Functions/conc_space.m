function [Qb1, Qb2] = conc_space(V1, V2)

dim1 = size(V1,2);
dim2 = size(V2,2);
dim = min([dim1, dim2]);
[U,~,~] = svd(full([V1, V2]));
V12 = U(:, 1:dim);
Qb1 = (V1'*V1)\(V1'*V12);
Qb2 = (V2'*V2)\(V2'*V12);