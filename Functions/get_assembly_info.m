function assembly_info = get_assembly_info(UC_dofs, pattern_dofs, pattern_uc)

nLR = numel(UC_dofs.L);
nBT = numel(UC_dofs.B);
nC = numel(UC_dofs.BL);
nI = numel(UC_dofs.I);
nMax = max([nLR,nBT,nC,nI]);
nSet = sum(pattern_dofs>0, 'all');
assembly_info.sn_data = zeros(nSet, nMax+1);

%% Map from sets to dofs
dof_id = 0;
for iSet = 1:size(pattern_dofs,1)
    for jSet = 1:size(pattern_dofs,2)
        if pattern_dofs(iSet,jSet) ~= 0
            set_id = pattern_dofs(iSet,jSet);
            switch mod(iSet,2) + mod(jSet,2)*1i
                case 0+0*1i
                    nNEW = nI;
                case 0+1i
                    nNEW = nLR;
                case 1+0*1i
                    nNEW = nBT;  
                case 1+1i
                    nNEW = nC; 
            end
            assembly_info.sn_data(set_id,1:nNEW+1) = [set_id, dof_id+1:dof_id+nNEW];
            dof_id = dof_id + nNEW;
        end
    end
end

%% Map from UCs to sets
nUC = sum(pattern_uc>0, 'all');
assembly_info.us_data = zeros(nUC, 10);
for iUC = 1:size(pattern_uc,1)
    for jUC = 1:size(pattern_uc,2)
        if pattern_uc(iUC,jUC) ~= 0
            uc_id = pattern_uc(iUC, jUC);
            iSet = 2*iUC-1:2*iUC+1;
            jSet = 2*jUC-1:2*jUC+1;
            lcset = pattern_dofs(iSet,jSet);
            assembly_info.us_data(uc_id,:) = [uc_id, lcset(2,2), lcset(2,1), lcset(2,3)...
                                                     lcset(3,2), lcset(1,2), lcset(3,1)...
                                                     lcset(3,3), lcset(1,3), lcset(1,1)];
        end
    end
end

assembly_info.reordering = [UC_dofs.I, UC_dofs.L, UC_dofs.R...
                            UC_dofs.B, UC_dofs.T, UC_dofs.BL...
                            UC_dofs.BR, UC_dofs.TR, UC_dofs.TL];
assembly_info.UC_nDOF = UC_dofs.nDOF;