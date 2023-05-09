function x_resp = EMD(K, M, F, R, freq, nModeI)

%% initialize the reduced system
if nModeI < size(K,1)/1.5
    disp("EMD...")
    [V, ~] = eigs(real(K+K')/2, real(M+M')/2, nModeI, 'smallestabs');
    stt_mode = real(K+K')\real(F);
    V = [V, stt_mode/norm(stt_mode)];
    K_red = V'*K*V;
    M_red = V'*M*V;
    F_red = V'*F;
    R_red = V'*R;
else
    K_red = K;
    M_red = M;
    F_red = F;
    R_red = R;
end

%% Frequency sweep
for iFreq = 1:numel(freq)
    if mod(iFreq,100) == 0
        fprintf('Freq: %d, ', iFreq);
    end
    omega = 2*pi*freq(iFreq);
    D_red = -omega^2*M_red + K_red;
    x_red = linsolve(D_red,F_red);
    x_resp(iFreq) = R_red'*x_red;
end

