function V_connected = compute_assemblied_basis(uc_array, assembly_info)

sn_data = assembly_info.sn_data;
us_data = assembly_info.us_data;
reordering = assembly_info.reordering;
nDOF_global = max(sn_data,[],'all');
UC_DOF = assembly_info.UC_nDOF;
nBasis = size(uc_array{1}.Basis, 2);
V_connected = zeros(nDOF_global,nBasis);

%% Generate assembled basis
krow = [1:UC_DOF];
for iUC = 1:size(us_data,1)
    sets_local = us_data(iUC,2:end);
    dofs_local = nonzeros(sn_data(sets_local,2:end)');
    V_connected(dofs_local, :) = uc_array{iUC}.Basis(reordering, :);
end


