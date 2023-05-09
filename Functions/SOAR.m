function [x_resp] = SOAR(K, M, F, R, nMode, freq, tol)

x_resp = [];
%% initialize the reduced system
V = [];
V = compute_basis2(M, K, F, V, 2*pi*500, 300);
% V = compute_basis2(M, K, F, V, 2*pi*800, min([ceil(nMode/3),60]));
% V = compute_basis2(M, K, F, V, 2*pi*900, min([ceil(nMode/3),40]));

[M_red, K_red, F_red] = compute_projection(M, K, F, V);
id_check = round(rand(100,1)*length(freq));
%% Frequency sweep
% error_array = [];
% error_array_ = [];
for iFreq = 1:numel(freq)
    if mod(iFreq,100) == 0
        fprintf('Freq: %d, ', iFreq);
    end
    omega = 2*pi*freq(iFreq);
    D_red = -omega^2*M_red + K_red;
    x_red = linsolve(D_red,F_red);
    x = V*x_red;
%     error_array = [error_array, norm((-omega^2*M + K)\F-x)/norm(x)];
%     error_array_ = [error_array_, norm((-omega^2*M + K)*x-F)/norm(F)];
    if sum(iFreq==id_check)>0
        error = norm((-omega^2*M + K)*x-F)/norm(F);
        if error > tol
            V = compute_basis2(M, K, F, V, omega*1.2, 300);
            [M_red, K_red, F_red] = compute_projection(M, K, F, V);
            D_red = -omega^2*M_red + K_red;
            x_red = linsolve(D_red,F_red);
            x = V*x_red;
        end
    end
    x_resp = [x_resp, R'*x];
end
% figure
% semilogy(freq, error_array)
% figure
% semilogy(freq, error_array_)
% a = 1;

function V = compute_basis(M, K, F, V, mMode)
dK = decomposition(K,'lu');
r = dK\F;
R0 = norm(r);
fprintf('\n');
for jMode = 1:mMode
    if mod(jMode,20) == 0
        fprintf('Mode: %d, ', jMode);
    end
    % add modes
    if jMode > 1
        tic
        r = dK\(M*r);
        toc
        tic
        r = linsolve(K,(M*r));
        toc
    end    
    R = norm(r);
    if R/R0 < 1e-12
        break;
    end
    r = r/R;
    V(:,end+1) = r;
%     V = orth(V);
    norm(V'*V-eye(size(V,2)))
end
fprintf('\n');
V = orth(V);

function V = compute_basis2(M, K, F, V, omega, mMode)
flag_sp = issparse(M);
A = real(K - omega^2*M);
B = real(- omega*2*M);
C = real(- M);
if flag_sp
    [L,U,P,Q] = lu(A);
else
    dA = decomposition(A,'lu');
end
p = 0*F;
if flag_sp
    q = Q*(U\(L\(P*F)));
else
    q = dA\F;
end
R0 = norm(q);
fprintf('\n');
for jMode = 1:mMode
    if mod(jMode,20) == 0
        fprintf('Mode: %d, ', jMode);
    end
    % add modes
    if jMode > 1
        if flag_sp
            r = Q*(U\(L\(P*(C*p + B*q))));
        else
            r = dA\(C*p + B*q);
        end
        for iMode = 1:size(V,2)
            r = r - V(:,iMode)* (V(:,iMode)'*r);
        end
        R = norm(r);
        if R/R0 < 1e-16
            break;
        end
        r = r/R;
        V(:,end+1) = r;
%         V = orth(V);
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
%         V = orth(V);
    end    
end
fprintf('\n');
% V = orth(V);

function [M_red, K_red, F_red] = compute_projection(M, K, F, V)
M_red = V'*M*V;
K_red = V'*K*V;
F_red = V'*F;

