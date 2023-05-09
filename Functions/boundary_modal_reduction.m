function [M_red,K_red,C_red,Dofs_red,V,L_A] = boundary_modal_reduction(M,K,C,Dofs,M_redI,K_redI,Dofs_redI,n_modes_A, Phi_I, Psi_IA, method_ort, tol_SVD)

%% Verify if n_modes_I is not larger than Dofs.I
if n_modes_A>length(Dofs_redI.A)
    error('n_modes_A is larger than amount of boundary dofs!')
end

%% Determine matrix partitions
M_AA = real(M_redI(Dofs_redI.A,Dofs_redI.A));
K_AA = real(K_redI(Dofs_redI.A,Dofs_redI.A));

% Determine sub-indices of A-portion; sorted [L R B T BL BR TR TL] already
% in interior modal reduction
Dofs_red_A.L = [1:length(Dofs_redI.L)];
Dofs_red_A.R = [Dofs_red_A.L(end)+1:Dofs_red_A.L(end)+length(Dofs_redI.R)];
Dofs_red_A.B = [Dofs_red_A.R(end)+1:Dofs_red_A.R(end)+length(Dofs_redI.B)];
Dofs_red_A.T = [Dofs_red_A.B(end)+1:Dofs_red_A.B(end)+length(Dofs_redI.T)];
Dofs_red_A.BL = [Dofs_red_A.T(end)+1:Dofs_red_A.T(end)+length(Dofs_redI.BL)];
Dofs_red_A.BR = [Dofs_red_A.BL(end)+1:Dofs_red_A.BL(end)+length(Dofs_redI.BR)];
Dofs_red_A.TR = [Dofs_red_A.BR(end)+1:Dofs_red_A.BR(end)+length(Dofs_redI.TR)];
Dofs_red_A.TL = [Dofs_red_A.TR(end)+1:Dofs_red_A.TR(end)+length(Dofs_redI.TL)];

%% Calculate normal mode shapes describing boundary motion
% Delete potential zero rows/columns (to avoid singularities)
delete_rows_A = ~any(K_AA,2);
K_AA(delete_rows_A,:) = [];
K_AA(:,delete_rows_A) = [];
M_AA(delete_rows_A,:) = [];
M_AA(:,delete_rows_A) = [];
keep_rows_A = setdiff(1:length(Dofs_redI.A),find(delete_rows_A));

% Calculate boundary normal modes
Phi_A = zeros(length(Dofs_redI.A),n_modes_A); %initialize
[Phi_A(keep_rows_A,:),~] = eigs(K_AA,M_AA,n_modes_A,0);
% [Phi_A,lambda_n_A] = eigs(K_AA,M_AA,n_modes_A,0);
% omega_n_A = (diag(lambda_n_A).^0.5)/(2*pi);

%% Calculate combined boundary mode sets LR, BT, BL-BR-TR-TL
% LR
% Combine mode shapes
Phi_LR = [Phi_A(Dofs_red_A.L,:),Phi_A(Dofs_red_A.R,:)];
% Remove zero rows LR (to avoid numerical noise from orthogonalization)
delete_rows_LR = ~any(Phi_LR,2);
keep_rows_LR = setdiff(1:length(Phi_LR(:,1)),find(delete_rows_LR));
% keep_rows_LR = 1:length(Phi_LR(:,1)); % keep all rows
% Verify if n_cols >= n_rows: if so --> replace by identity matrix;
% otherwise carry on with Phi_LR and orthogonalize.
if length(Phi_LR(1,:))>=length(Phi_LR(:,1))
    Phi_LR = eye(length(Phi_LR(:,1)));
else % orthogonalize modeset
    switch method_ort
        case 'SVD'
            [U_LR,S_LR,~] = svd(Phi_LR(keep_rows_LR,:));
            % Phi_LR_orth = U_LR(:,1:length(Phi_LR(1,:)));
            % Retain non-null singular values (defined by tol_SVD)
            % according to range-space tolerance (w.r.t. max singular
            % value)
            Phi_LR_orth = U_LR(:,(diag(S_LR)/S_LR(1,1))>tol_SVD);
        case 'QR'
            [Q_LR,~,~] = qr(Phi_LR(keep_rows_LR,:));
            Phi_LR_orth = Q_LR(:,1:length(Phi_LR(1,:)));
        otherwise
            error('Invalid orthogonalization method')
    end
    Phi_LR = zeros(length(Phi_LR(:,1)),length(Phi_LR_orth(1,:)));
    Phi_LR(keep_rows_LR,:) = Phi_LR_orth;
end

% BT
% Combine mode shapes
Phi_BT = [Phi_A(Dofs_red_A.B,:),Phi_A(Dofs_red_A.T,:)];
% Remove zero rows BT (to avoid numerical noise from orthogonalization)
delete_rows_BT = ~any(Phi_BT,2);
keep_rows_BT = setdiff(1:length(Phi_BT(:,1)),find(delete_rows_BT));
% keep_rows_BT = 1:length(Phi_BT(:,1)); % keep all rows
% Verify if n_cols >= n_rows: if so --> replace by identity matrix;
% otherwise carry on with Phi_LR and orthogonalize.
if length(Phi_BT(1,:))>=length(Phi_BT(:,1))
    Phi_BT = eye(length(Phi_BT(:,1)));
else % orthogonalize modeset
    switch method_ort
        case 'SVD'
            [U_BT,S_BT,~] = svd(Phi_BT(keep_rows_BT,:));
            % Phi_BT_orth = U_BT (:,1:length(Phi_BT(1,:)));
            % Retain non-null singular values (defined by tol_SVD)
            % according to range-space tolerance (w.r.t. max singular
            % value)
            Phi_BT_orth = U_BT(:,(diag(S_BT)/S_BT(1,1))>tol_SVD);
        case 'QR'
            [Q_BT,~] = qr(Phi_BT(keep_rows_BT,:));
            Phi_BT_orth = Q_BT(:,1:length(Phi_BT(1,:)));
        otherwise
            error('Invalid orthogonalization method')
    end
    Phi_BT = zeros(length(Phi_BT(:,1)),length(Phi_BT_orth(1,:)));
    Phi_BT(keep_rows_BT,:) = Phi_BT_orth;
end

% BL-BR-TR-TL
% Combine mode shapes
Phi_C = [Phi_A(Dofs_red_A.BL,:),Phi_A(Dofs_red_A.BR,:),Phi_A(Dofs_red_A.TR,:),Phi_A(Dofs_red_A.TL,:)];
% Remove zero rows C (to avoid numerical noise from orthogonalization)
delete_rows_C = ~any(Phi_C,2);
keep_rows_C = setdiff(1:length(Phi_C(:,1)),find(delete_rows_C));
% keep_rows_C = 1:length(Phi_C(:,1)); % keep all rows
% Verify if n_cols >= n_rows: if so --> replace by identity matrix;
% otherwise carry on with Phi_LR and orthogonalize.
if length(Phi_C(1,:))>=length(Phi_C(:,1))
    Phi_C = eye(length(Phi_C(:,1)));
else % orthogonalize modeset
    switch method_ort
        case 'SVD'
            [U_C,S_C,~] = svd(Phi_C(keep_rows_C,:));
            % Phi_C_orth = U_C(:,1:length(Phi_C(1,:)));
            % Retain non-null singular values (defined by tol_SVD)
            % according to range-space tolerance (w.r.t. max singular
            % value)
            Phi_C_orth = U_C(:,(diag(S_C)/S_C(1,1))>tol_SVD);
        case 'QR'
            [Q_C,~] = qr(Phi_C(keep_rows_C,:));
            Phi_C_orth = Q_C(:,1:length(Phi_C(1,:)));
        otherwise
            error('Invalid orthogonalization method')
    end
    Phi_C = zeros(length(Phi_C(:,1)),length(Phi_C_orth(1,:)));
    Phi_C(keep_rows_C,:) = Phi_C_orth;
end

%% Create boundary transformation matrix
% % L = sparse(2*length(Phi_LR(:,1))+2*length(Phi_BT(:,1))+4*length(Phi_C(:,1)),2*length(Phi_LR(1,:))+2*length(Phi_BT(1,:))+4*length(Phi_C(1,:)),0);
L_A = zeros(2*length(Phi_LR(:,1))+2*length(Phi_BT(:,1))+4*length(Phi_C(:,1)),2*length(Phi_LR(1,:))+2*length(Phi_BT(1,:))+4*length(Phi_C(1,:)));
%LR
L_A(Dofs_red_A.L,1:length(Phi_LR(1,:))) = Phi_LR;
L_A(Dofs_red_A.R,length(Phi_LR(1,:))+1:2*length(Phi_LR(1,:))) = Phi_LR;
%BT
L_A(Dofs_red_A.B,2*length(Phi_LR(1,:))+1:2*length(Phi_LR(1,:))+length(Phi_BT(1,:))) = Phi_BT;
L_A(Dofs_red_A.T,2*length(Phi_LR(1,:))+length(Phi_BT(1,:))+1:2*length(Phi_LR(1,:))+2*length(Phi_BT(1,:))) = Phi_BT;
%BL-BR-TR-TL
L_A(Dofs_red_A.BL,2*length(Phi_LR(1,:))+2*length(Phi_BT(1,:))+1:2*length(Phi_LR(1,:))+2*length(Phi_BT(1,:))+length(Phi_C(1,:))) = Phi_C;
L_A(Dofs_red_A.BR,2*length(Phi_LR(1,:))+2*length(Phi_BT(1,:))+length(Phi_C(1,:))+1:2*length(Phi_LR(1,:))+2*length(Phi_BT(1,:))+2*length(Phi_C(1,:))) = Phi_C;
L_A(Dofs_red_A.TR,2*length(Phi_LR(1,:))+2*length(Phi_BT(1,:))+2*length(Phi_C(1,:))+1:2*length(Phi_LR(1,:))+2*length(Phi_BT(1,:))+3*length(Phi_C(1,:))) = Phi_C;
L_A(Dofs_red_A.TL,2*length(Phi_LR(1,:))+2*length(Phi_BT(1,:))+3*length(Phi_C(1,:))+1:2*length(Phi_LR(1,:))+2*length(Phi_BT(1,:))+4*length(Phi_C(1,:))) = Phi_C;

%% Calculate new reduced dofs
Dofs_red.I = Dofs_redI.I;
Dofs_red.L = [Dofs_red.I(end)+1:Dofs_red.I(end)+length(Phi_LR(1,:))];
Dofs_red.R = [Dofs_red.L(end)+1:Dofs_red.L(end)+length(Phi_LR(1,:))];
Dofs_red.B = [Dofs_red.R(end)+1:Dofs_red.R(end)+length(Phi_BT(1,:))];
Dofs_red.T = [Dofs_red.B(end)+1:Dofs_red.B(end)+length(Phi_BT(1,:))];
Dofs_red.BL = [Dofs_red.T(end)+1:Dofs_red.T(end)+length(Phi_C(1,:))];
Dofs_red.BR = [Dofs_red.BL(end)+1:Dofs_red.BL(end)+length(Phi_C(1,:))];
Dofs_red.TR = [Dofs_red.BR(end)+1:Dofs_red.BR(end)+length(Phi_C(1,:))];
Dofs_red.TL = [Dofs_red.TR(end)+1:Dofs_red.TR(end)+length(Phi_C(1,:))];
Dofs_red.A = [Dofs_red.L,Dofs_red.R,Dofs_red.B,Dofs_red.T,Dofs_red.BL,Dofs_red.BR,Dofs_red.TR,Dofs_red.TL];

%% Modify BMS transformation matrix to incorporate L_A
% Initialize
V = zeros(length(Dofs.I)+length(Dofs.A),length(Dofs_red.I)+length(Dofs_red.A));
% Assign to corresponding DOFs
V(Dofs.I,Dofs_red.I) = Phi_I;
V(Dofs.I,Dofs_red.A) = Psi_IA*L_A;
V(Dofs.A,Dofs_red.A) =  L_A;

%% Reduce original UC matrices
M_red = V.'*M*V;
K_red = V.'*K*V;
C_red = V.'*C*V;

%% Symmetrize
M_red = (M_red+M_red.')/2;
K_red = (K_red+K_red.')/2;
C_red = (C_red+C_red.')/2;

end