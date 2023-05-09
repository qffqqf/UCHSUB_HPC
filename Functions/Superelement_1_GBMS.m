classdef Superelement_1_GBMS

    properties
        F_UC_red
        R_UC_red
        M_UC_red
        K_UC_red
        U_UC_red
        V_UC
        V_UC_fom
        Basis_red
        Basis
        F_UC
        R_UC
        M_UC
        K_UC
        nModeI
        nModeA
        dofs
        dofs_fom
        dofs_red
        Proj
        Proj_delete
    end

    methods

        function obj = init_comp(obj)
            obj.Proj = speye(obj.dofs_red.nDOF);
            obj.Proj_delete = [];
        end

        function obj = update_boundary(obj, Q_b, side)
            nVb = size(Q_b, 2);
            switch side
                case "B"
                    obj.Proj(obj.dofs_red.B, obj.dofs_red.B(1):obj.dofs_red.B(1)-1+nVb) = Q_b;
                    obj.Proj_delete = [obj.Proj_delete, obj.dofs_red.B(1)+nVb:obj.dofs_red.B(end)];
                case "R"
                    obj.Proj(obj.dofs_red.R, obj.dofs_red.R(1):obj.dofs_red.R(1)-1+nVb) = Q_b;
                    obj.Proj_delete = [obj.Proj_delete, obj.dofs_red.R(1)+nVb:obj.dofs_red.R(end)];
                case "T"
                    obj.Proj(obj.dofs_red.T, obj.dofs_red.T(1):obj.dofs_red.T(1)-1+nVb) = Q_b;
                    obj.Proj_delete = [obj.Proj_delete, obj.dofs_red.T(1)+nVb:obj.dofs_red.T(end)];
                case "L"
                    obj.Proj(obj.dofs_red.L, obj.dofs_red.L(1):obj.dofs_red.L(1)-1+nVb) = Q_b;
                    obj.Proj_delete = [obj.Proj_delete, obj.dofs_red.L(1)+nVb:obj.dofs_red.L(end)];
            end
        end

        function obj = Superelement_1_GBMS(M_UC, K_UC, F_UC, R_UC, dofs)
            obj.M_UC = M_UC;
            obj.K_UC = K_UC;
            obj.F_UC = F_UC;
            obj.R_UC = R_UC;
            obj.dofs = dofs;
            obj.dofs.nDOF_A = numel(dofs.A);
            obj.dofs_fom = obj.dofs;
        end
         
        function obj = compute_modes(obj, nModeI, nModeA, obj_temp, freq)
            obj.nModeI = nModeI;
            obj.nModeA = nModeA;
            if isempty(obj_temp)
                if nModeI > numel(obj.dofs.I)-1 % rom>fom
                    obj.V_UC = speye(obj.dofs.nDOF);
                    obj.V_UC_fom = obj.V_UC(obj.dofs.A, :);
                    obj.dofs_red = obj.dofs;
                else
                    fprintf("Boundary modal reduction:1...")
                    K_II = obj.K_UC(obj.dofs.I,obj.dofs.I);
                    K_IA = obj.K_UC(obj.dofs.I,obj.dofs.A);
                    M_II = obj.M_UC(obj.dofs.I,obj.dofs.I);
                    M_II = real(M_II + M_II')/2;
                    K_II = real(K_II + K_II')/2;
                    stat_modes = -linsolve(full(K_II), full(K_IA));
                    transf = [eye(obj.dofs.nDOF_A);stat_modes];
                    fprintf('%d/%d\n', nModeA, obj.dofs.nDOF_A);
                    %% get boundary modes
                    M_AA = transf'*full(obj.M_UC)*transf;
                    K_AA = transf'*full(obj.K_UC)*transf;
                    M_AA = real(M_AA + M_AA')/2;
                    K_AA = real(K_AA + K_AA')/2;
                    [Phi_A,~] = eigs(K_AA,M_AA,nModeA,'smallestabs');
                    Phi_L = orth([Phi_A(obj.dofs.L,:),Phi_A(obj.dofs.R,:)]);
                    if size(Phi_L,2) > size(Phi_L,1)/1.2
                        Phi_L = eye(size(Phi_L,1));
                    end
                    nL = size(Phi_L,2);
                    Phi_R = orth([Phi_A(obj.dofs.L,:),Phi_A(obj.dofs.R,:)]);
                    if size(Phi_R,2) > size(Phi_R,1)/1.2
                        Phi_R = eye(size(Phi_R,1));
                    end
                    nR = size(Phi_R,2);
                    Phi_B = orth([Phi_A(obj.dofs.B,:),Phi_A(obj.dofs.T,:)]);
                    if size(Phi_B,2) > size(Phi_B,1)/1.2
                        Phi_B = eye(size(Phi_B,1));
                    end
                    nB = size(Phi_B,2);
                    Phi_T = orth([Phi_A(obj.dofs.B,:),Phi_A(obj.dofs.T,:)]);
                    if size(Phi_T,2) > size(Phi_T,1)/1.2
                        Phi_T = eye(size(Phi_T,1));
                    end
                    nT = size(Phi_T,2);
                    Phi_C = eye(numel(obj.dofs.BL));
                    nC = numel(obj.dofs.BL);
                    %% get reduced dofs index
                    obj.dofs_red.BL = [1:nC]; 
                    obj.dofs_red.B = obj.dofs_red.BL(end) + [1:nB];
                    obj.dofs_red.BR = obj.dofs_red.B(end) + [1:nC];
                    obj.dofs_red.R = obj.dofs_red.BR(end) + [1:nR];
                    obj.dofs_red.TR = obj.dofs_red.R(end) + [1:nC];
                    obj.dofs_red.T = obj.dofs_red.TR(end) + [1:nT];
                    obj.dofs_red.TL = obj.dofs_red.T(end) + [1:nC];
                    obj.dofs_red.L = obj.dofs_red.TL(end) + [1:nL];
                    obj.dofs_red.I = obj.dofs_red.L(end) + [1:nModeI];
                    obj.dofs_red.A = [obj.dofs_red.BL, obj.dofs_red.B, obj.dofs_red.BR, obj.dofs_red.R,...
                                      obj.dofs_red.TR, obj.dofs_red.T, obj.dofs_red.TL, obj.dofs_red.L];
                    obj.dofs_red.nDOF_A = nL + nR + nB + nT + nC*4;
                    obj.dofs_red.nDOF = obj.dofs_red.nDOF_A + nModeI;
                    %% get split modes
                    L_A = zeros(obj.dofs.nDOF_A, obj.dofs_red.nDOF_A);
                    L_A(obj.dofs.BL, obj.dofs_red.BL) = Phi_C;
                    L_A(obj.dofs.B, obj.dofs_red.B) = Phi_B;
                    L_A(obj.dofs.BR, obj.dofs_red.BR) = Phi_C;
                    L_A(obj.dofs.R, obj.dofs_red.R) = Phi_R;
                    L_A(obj.dofs.TR, obj.dofs_red.TR) = Phi_C;
                    L_A(obj.dofs.T, obj.dofs_red.T) = Phi_T;
                    L_A(obj.dofs.TL, obj.dofs_red.TL) = Phi_C;
                    L_A(obj.dofs.L, obj.dofs_red.L) = Phi_L;
                    fprintf("Interior modal reduction:1...")
                    %% get interior modes
                    fprintf('%d/%d\n', nModeI, numel(obj.dofs.I));
                    %int_modes = krylov(M_II, K_II, ones(size(K_II,1),1), [], freq, nModeI);
                    [int_modes,~] = eigs(K_II, M_II, nModeI, 'smallestabs');
                    obj.nModeI = nModeI;
                    %% Get subspace
                    obj.V_UC = zeros(obj.dofs.nDOF, obj.dofs_red.nDOF);
                    obj.V_UC(obj.dofs.A, obj.dofs_red.A) = L_A;
                    obj.V_UC(obj.dofs.I, obj.dofs_red.A) = stat_modes*L_A;
                    obj.V_UC(obj.dofs.I, obj.dofs_red.I) = int_modes;
                    obj.V_UC_fom = obj.V_UC(obj.dofs.A, :);
                end
            else
                obj.V_UC = obj_temp.V_UC;
                obj.nModeI = obj_temp.nModeI;
                obj.nModeA = obj_temp.nModeA;
                obj.dofs_red = obj_temp.dofs_red;
                obj.V_UC_fom = obj_temp.V_UC_fom;
            end
        end

        function obj = projection(obj, obj_temp)
            if isempty(obj_temp)
                disp("projection:1...")
                obj.M_UC_red = obj.V_UC'* full(obj.M_UC)* obj.V_UC;
                obj.K_UC_red = obj.V_UC'* full(obj.K_UC)* obj.V_UC;
                obj.M_UC_red = real(full(obj.M_UC_red + obj.M_UC_red')/2);
                obj.K_UC_red = real(full(obj.K_UC_red + obj.K_UC_red')/2);
                obj.F_UC_red = obj.V_UC'* obj.F_UC;
                obj.R_UC_red = obj.V_UC'* obj.R_UC;
            else
                obj.M_UC_red = obj_temp.M_UC_red;
                obj.K_UC_red = obj_temp.K_UC_red;
                obj.F_UC_red = obj_temp.V_UC'* obj.F_UC;
                obj.R_UC_red = obj_temp.V_UC'* obj.R_UC;
            end
        end

        function obj = free_memory(obj)
            obj.M_UC_red = [];
            obj.K_UC_red = [];
            obj.F_UC_red = [];
            obj.R_UC_red = [];
            
            obj.M_UC = [];
            obj.K_UC = [];
            obj.F_UC = [];
            obj.R_UC = [];
            obj.V_UC = [];
            obj.V_UC_fom = [];
        end

        function obj = recovery(obj)
            obj.Basis = obj.V_UC* obj.Basis_red;
        end

    end
end