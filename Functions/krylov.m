function V = krylov(M, K, F, V, omega, mMode)
A = real(K - omega^2*M);
B = real(- omega*2*M);
C = real(- M);
% A = (K - omega^2*M);
% B = (- omega*2*M);
% C = (- M);
dA = decomposition(A,'lu');
p = 0*F;
q = dA\F;
R0 = norm(q);
for jMode = 1:mMode
    % add modes
    if jMode > 1
        r = dA\(C*p + B*q);
        for iMode = 1:size(V,2)
            r = r - V(:,iMode)* (V(:,iMode)'*r);
        end
        R = norm(r);
        if R/R0 < 1e-16
            break;
        end
        r = r/R;
        V(:,end+1) = r;
        p = q;
        q = r;
    else
        r = q;
        for iMode = 1:size(V,2)
            r = r - V(:,iMode)* (V(:,iMode)'*r);
        end
        R = norm(r);
        if R/R0 < 1e-16
            break;
        end
        r = r/R;
        V(:,end+1) = r;
    end    
end