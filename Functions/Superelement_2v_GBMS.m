classdef Superelement_2v_GBMS

    properties
        F_UC_red
        R_UC_red
        M_UC_red
        K_UC_red
        U_UC_red
        V_UC
        Basis_red
        Basis
        F_UC
        R_UC
        M_UC
        K_UC
        nModeI
        nModeA
        dofs
        dofs_red
        ev_map
    end

    methods

        function [obj, children] = find_map(obj, ele_pack, obj_temp)
            if isempty(obj_temp)
                %% init
                ele_bl = ele_pack{1};
                ele_tl = ele_pack{2};
                %% align dofs
                dofs_BL = [1:numel(ele_bl.dofs_red.BL)];
                dofs_B = dofs_BL(end) + [1:numel(ele_bl.dofs_red.B)];
                dofs_BR = dofs_B(end) + [1:numel(ele_bl.dofs_red.BR)];
                dofs_R1 = dofs_BR(end) + [1:numel(ele_bl.dofs_red.R)];
                dofs_R2 = dofs_R1(end) + [1:numel(ele_bl.dofs_red.BR)];
                dofs_R3 = dofs_R2(end) + [1:numel(ele_tl.dofs_red.R)];
                dofs_TR = dofs_R3(end) + [1:numel(ele_tl.dofs_red.TR)];
                dofs_T = dofs_TR(end) + [1:numel(ele_tl.dofs_red.T)];
                dofs_TL = dofs_T(end) + [1:numel(ele_tl.dofs_red.TL)];
                dofs_L1 = dofs_TL(end) + [1:numel(ele_bl.dofs_red.L)];
                dofs_L2 = dofs_L1(end) + [1:numel(ele_bl.dofs_red.TL)];
                dofs_L3 = dofs_L2(end) + [1:numel(ele_tl.dofs_red.L)];
                nI1 = numel(ele_tl.dofs_red.B);
                dofs_I1 = dofs_L3(end) + [1:nI1];
                dofs_I2 = dofs_I1(end) + [1:numel(ele_bl.dofs_red.I)];
                dofs_I3 = dofs_I2(end) + [1:numel(ele_tl.dofs_red.I)];
                obj.dofs.BL = dofs_BL;
                obj.dofs.B = dofs_B;
                obj.dofs.BR = dofs_BR;
                obj.dofs.R = [dofs_R1,dofs_R2,dofs_R3];
                obj.dofs.TR = dofs_TR;
                obj.dofs.T = dofs_T;
                obj.dofs.TL = dofs_TL;
                obj.dofs.L = [dofs_L1,dofs_L2,dofs_L3];
                obj.dofs.I = [dofs_I1,dofs_I2,dofs_I3];
                obj.dofs.A = [obj.dofs.BL,obj.dofs.B,obj.dofs.BR,obj.dofs.R,obj.dofs.TR,obj.dofs.T,obj.dofs.TL,obj.dofs.L];
                obj.dofs.nDOF = numel([obj.dofs.A, obj.dofs.I]);
                obj.dofs.nDOF_A = numel(obj.dofs.A);
                obj.ev_map = cell(2,1);
                obj.ev_map{1} = [dofs_BL,dofs_B,dofs_BR,dofs_R1,dofs_R2,dofs_I1,dofs_L2,dofs_L1,dofs_I2];
                obj.ev_map{2} = [dofs_L2,dofs_I1,dofs_R2,dofs_R3,dofs_TR,dofs_T,dofs_TL,dofs_L3,dofs_I3];
                children = {ele_bl, ele_tl};
            else
                ele_bl = ele_pack{1};
                ele_tl = ele_pack{2};
                obj.dofs = obj_temp.dofs;
                obj.ev_map = obj_temp.ev_map;
                children = {ele_bl, ele_tl};
            end
        end
         
        function obj = assembly(obj, ele_pack, obj_temp)
            if isempty(obj_temp)
                disp("Assembly:2...")
                obj.M_UC = zeros(obj.dofs.nDOF, obj.dofs.nDOF);
                obj.K_UC = zeros(obj.dofs.nDOF, obj.dofs.nDOF);
                for iEle = 1:2
                    dofs_local = nonzeros(obj.ev_map{iEle});
                    Isps = dofs_local;
                    obj.M_UC(Isps,Isps) = obj.M_UC(Isps,Isps) + ele_pack{iEle}.M_UC_red;
                    obj.K_UC(Isps,Isps) = obj.K_UC(Isps,Isps) + ele_pack{iEle}.K_UC_red;
                end
%                 obj.M_UC = real(obj.M_UC+obj.M_UC')/2;
%                 obj.K_UC = real(obj.K_UC+obj.K_UC')/2;
            else
                obj.M_UC = obj_temp.M_UC;
                obj.K_UC = obj_temp.K_UC;
            end

            obj.F_UC = zeros(obj.dofs.nDOF, 1);
            obj.R_UC = zeros(obj.dofs.nDOF, 1);
            for iEle = 1:2
                dofs_local = nonzeros(obj.ev_map{iEle});
                Isps = dofs_local;
                obj.F_UC(Isps) = obj.F_UC(Isps) + ele_pack{iEle}.F_UC_red;
                obj.R_UC(Isps) = obj.R_UC(Isps) + ele_pack{iEle}.R_UC_red;
            end
            obj.nModeI = round(mean([ele_pack{1}.nModeI, ele_pack{2}.nModeI])*1.7); %% change rom here!
        end

        function obj = compute_modes(obj, obj_temp, freq)
            if isempty(obj_temp)
                fprintf('%d/%d/%d\n', numel(obj.dofs.I), numel(obj.dofs.A), obj.nModeI);
                if or(numel(obj.dofs.I) < numel(obj.dofs.A), numel(obj.dofs.I) < obj.nModeI*1.5)
                    obj.V_UC = speye(obj.dofs.nDOF);
                    obj.dofs_red = obj.dofs;
                else
                    K_II = obj.K_UC(obj.dofs.I,obj.dofs.I);
                    K_IA = obj.K_UC(obj.dofs.I,obj.dofs.A);
                    M_II = obj.M_UC(obj.dofs.I,obj.dofs.I);
%                     M_II = real(M_II + M_II')/2;
%                     K_II = real(K_II + K_II')/2;
                    stat_modes = -linsolve(full(K_II), full(K_IA));
                    %% get reduced dofs index
                    obj.dofs_red = obj.dofs;
                    obj.dofs_red.I = obj.dofs_red.L(end) + [1:obj.nModeI];
                    obj.dofs_red.nDOF = obj.dofs_red.nDOF_A + obj.nModeI;
                    %% get interior modes
                    fprintf("Interior modal reduction:4...")
                    fprintf('%d/%d\n', obj.nModeI, numel(obj.dofs.I));   
                    [int_modes,~] = eigs(K_II, M_II, obj.nModeI, 'smallestabs');
                    %% Get subspace
                    obj.V_UC = zeros(obj.dofs.nDOF, obj.dofs_red.nDOF);
                    obj.V_UC(obj.dofs.A, obj.dofs_red.A) = eye(obj.dofs_red.nDOF_A);
                    obj.V_UC(obj.dofs.I, obj.dofs_red.A) = stat_modes;
                    obj.V_UC(obj.dofs.I, obj.dofs_red.I) = int_modes;
                end
            else
                obj.dofs_red = obj_temp.dofs_red;
                obj.V_UC = obj_temp.V_UC;
                obj.nModeI = obj_temp.nModeI;
            end
        end

        function obj = projection(obj, obj_temp)
            if isempty(obj_temp)
                disp("Projection:4...")
                obj.M_UC_red = full(obj.V_UC'* obj.M_UC* obj.V_UC);
                obj.K_UC_red = full(obj.V_UC'* obj.K_UC* obj.V_UC);
%                 obj.M_UC_red = real(full(obj.M_UC_red + obj.M_UC_red')/2);
%                 obj.K_UC_red = real(full(obj.K_UC_red + obj.K_UC_red')/2);
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
        end

        function obj = recovery(obj)
            obj.Basis = obj.V_UC* obj.Basis_red;
        end

        function ele_pack = pass_info(obj, ele_pack)
            for iEle = 1:2
                dofs_local = nonzeros(obj.ev_map{iEle});
                ele_pack{iEle}.Basis_red = ele_pack{iEle}.Proj* obj.Basis(dofs_local,:);
            end
        end

    end
end