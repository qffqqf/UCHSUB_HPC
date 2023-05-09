function N = null_space(A, nMode)

[~,~,V] = svd(A);
N = null(A);
dim = size(N,2);
if dim > nMode
    N = N(:,1:nMode);
else
    count = 0;
    while size(N,2) ~= nMode
        q = V(:,end-count);
        count = count + 1;
        q = q - N* (N'*q);
        q = q/norm(q + eps);
        N(:,end+1) = q;
    end
end