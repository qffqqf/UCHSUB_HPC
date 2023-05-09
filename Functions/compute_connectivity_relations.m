function pattern_dofs = compute_connectivity_relations(pattern_uc)

%% Augmente pattern matrix
N_UC_x = size(pattern_uc, 1);
N_UC_y = size(pattern_uc, 2);
pattern_dofs = zeros(2*N_UC_x+1, 2*N_UC_y+1);
pattern_dofs(2:2:end,2:2:end) = 100*pattern_uc;

%% generate connections
for iUCx = 2:2:2*N_UC_x+1
    for iUCy = 2:2:2*N_UC_y+1
        if pattern_dofs(iUCx,iUCy) > 0
            pattern_dofs(iUCx-1:iUCx+1,iUCy-1:iUCy+1) = ...
                pattern_dofs(iUCx-1:iUCx+1,iUCy-1:iUCy+1) - ones(3,3);
        end
    end
end
uc_id = find(pattern_dofs>0);
interface_id = find(pattern_dofs<0);
pattern_dofs(uc_id) = 1:numel(uc_id);
pattern_dofs(interface_id) = numel(uc_id)+1:numel(uc_id)+numel(interface_id);

            