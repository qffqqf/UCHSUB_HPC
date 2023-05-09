function eig_freq = hgbms_init(uc_index, root_index, ...
                              F_UC, M_UC, K_UC, UC_dofs, pointer_mx, brothers_mx, nModeI, nModeA, nEig, freq)

% get UCs
disp("get uc")
for iIndex = numel(root_index):-1:1
    
    se_id = root_index(iIndex);
    uc_id = uc_index(iIndex);
    F_UC_ = zeros(size(F_UC));
    se_array{se_id} = Superelement_1_GBMS(M_UC, K_UC, F_UC_, UC_dofs);
    if sum(brothers_mx(se_id,:)) > 0
        fprintf('# %d; ',iIndex);
        brother_id = find(brothers_mx(se_id,:) == 1);
        se_array{se_id} = se_array{se_id}.compute_modes(nModeI, nModeA, se_array{brother_id}, freq);
        se_array{se_id} = se_array{se_id}.projection(se_array{brother_id});
    else
        fprintf('\n## UC: %d\n',iIndex);
        se_array{se_id} = se_array{se_id}.compute_modes(nModeI, nModeA, [], freq);
        se_array{se_id} = se_array{se_id}.projection([]);
    end
end

% Condensation
fprintf("\nCondensation\n")
pointer_mx_sum = sum(pointer_mx,1);
for iFather = numel(pointer_mx_sum):-1:1
    if pointer_mx_sum(iFather) == 4
        
        se_array{iFather} = Superelement_4_GBMS;
        children_index = find(pointer_mx(:,iFather));
        children = se_array(children_index);
        if sum(brothers_mx(iFather,:)) > 0
            fprintf('# %d; ',iFather);
            brother_id = find(brothers_mx(iFather,:) == 1);
            se_array{iFather} = se_array{iFather}.find_fom_map(children, se_array{brother_id});
            [se_array{iFather}, ~] = se_array{iFather}.find_map(children, se_array{brother_id});
            se_array{iFather} = se_array{iFather}.assembly(children, se_array{brother_id});
            se_array{iFather} = se_array{iFather}.compute_modes(se_array{brother_id}, freq);
            se_array{iFather} = se_array{iFather}.projection(se_array{brother_id});
        else
            if iFather > 1
                fprintf('\n#### Super UC: %d\n',iFather);
                se_array{iFather} = se_array{iFather}.find_fom_map(children, []);
                [se_array{iFather}, children] = se_array{iFather}.find_map(children, []);
                se_array{iFather} = se_array{iFather}.assembly(children, []);
                se_array{iFather} = se_array{iFather}.compute_modes([], freq);
                se_array{iFather} = se_array{iFather}.projection([]);
            else
                fprintf('\n#### Super UC: %d\n',iFather);
                se_array{iFather} = se_array{iFather}.find_fom_map(children, []);
                [se_array{iFather}, children] = se_array{iFather}.find_map(children, []);
                se_array{iFather} = se_array{iFather}.assembly(children, []);
            end
        end
    elseif pointer_mx_sum(iFather) == 2
        children_index = reshape(find(pointer_mx(:,iFather)),[2,1]);
        uc_id = uc_index(find(sum(root_index==children_index,1)));
        if abs(uc_id(1)-uc_id(2)) == 1
            se_array{iFather} = Superelement_2v_GBMS;
        else 
            se_array{iFather} = Superelement_2h_GBMS;
        end
        children = se_array(children_index);
        if sum(brothers_mx(iFather,:)) > 0
            fprintf('# %d; ',iFather);
            brother_id = find(brothers_mx(iFather,:) == 1);
            se_array{iFather} = se_array{iFather}.find_fom_map(children, se_array{brother_id});
            [se_array{iFather},~] = se_array{iFather}.find_map(children, se_array{brother_id});
            se_array{iFather} = se_array{iFather}.assembly(children, se_array{brother_id});
            se_array{iFather} = se_array{iFather}.compute_modes(se_array{brother_id}, freq);
            se_array{iFather} = se_array{iFather}.projection(se_array{brother_id});
        else
            fprintf('\n#### Super UC: %d\n',iFather);
            se_array{iFather} = se_array{iFather}.find_fom_map(children, []);
            [se_array{iFather},children] = se_array{iFather}.find_map(children, []);
            se_array{iFather} = se_array{iFather}.assembly(children, []);
            se_array{iFather} = se_array{iFather}.compute_modes([], freq);
            se_array{iFather} = se_array{iFather}.projection([]);
        end
    end
end
% solve
disp('Solving Eig ...')
tic
se_array{1} = se_array{1}.solve_eig(nEig);
toc
eig_freq = se_array{1}.eigen_freq;