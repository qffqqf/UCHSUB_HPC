function frf = hss(uc_model, param, tree_model, options)

%% UC reduction
disp("UC MOR...")
nSE = numel(tree_model.root_index); % SE: superelement
for iIndex = nSE:-1:1
    % Unfolding
    se_id = tree_model.root_index(iIndex);
    uc_id = tree_model.uc_index(iIndex);
    M_UC = uc_model.M_UC;
    K_UC = uc_model.K_UC;
    F_UC = param.forc(uc_id)* uc_model.F_UC;
    R_UC = param.resp(uc_id)* uc_model.R_UC;
    % Reduction
    se_array{se_id} = Superelement_1_GBMS(M_UC, K_UC, F_UC, R_UC, uc_model.UC_dofs);
    fprintf('\n# UC: %d\n',iIndex);
    if sum(tree_model.brothers_mx(se_id,:)) > 0
        brother_id = find(tree_model.brothers_mx(se_id,:) == 1);
        se_array{se_id} = se_array{se_id}.compute_modes(options.nModeI, options.nModeA, se_array{brother_id}, options.exppt);
        se_array{se_id} = se_array{se_id}.projection(se_array{brother_id});
    else
        se_array{se_id} = se_array{se_id}.compute_modes(options.nModeI, options.nModeA, [], options.exppt);
        se_array{se_id} = se_array{se_id}.projection([]);
    end

end

%% Multilevel SE reduction
fprintf("\nMultilevel MOR...\n")
pointers = sum(tree_model.pointer_mx, 1);
nPT = numel(pointers);
for iFather = nPT:-1:1
    if pointers(iFather) == 4
        se_array{iFather} = Superelement_4_GBMS;
        children_index = find(tree_model.pointer_mx(:,iFather));
        children = se_array(children_index);
        fprintf('\n## Super UC: %d\n',iFather);
        if sum(tree_model.brothers_mx(iFather,:)) > 0
            brother_id = find(tree_model.brothers_mx(iFather,:) == 1);
            [se_array{iFather}, children] = se_array{iFather}.find_map(children, se_array{brother_id});
            se_array{iFather} = se_array{iFather}.assembly(children, se_array{brother_id});
            se_array{iFather} = se_array{iFather}.compute_modes(se_array{brother_id}, options.exppt);
            se_array{iFather} = se_array{iFather}.projection(se_array{brother_id});
        else
            if iFather > 1
                [se_array{iFather}, children] = se_array{iFather}.find_map(children, []);
                se_array{iFather} = se_array{iFather}.assembly(children, []);
                se_array{iFather} = se_array{iFather}.compute_modes([], options.exppt);
                se_array{iFather} = se_array{iFather}.projection([]);
            else
                [se_array{iFather}, children] = se_array{iFather}.find_map(children, []);
                se_array{iFather} = se_array{iFather}.assembly(children, []);
            end
        end
        se_array(children_index) = children;
    elseif pointers(iFather) == 2
        children_index = reshape(find(tree_model.pointer_mx(:,iFather)),[2,1]);
        uc_id = tree_model.uc_index(find(sum(tree_model.root_index == children_index,1)));
        if abs(uc_id(1)-uc_id(2)) == 1
            se_array{iFather} = Superelement_2v_GBMS;
        else 
            se_array{iFather} = Superelement_2h_GBMS;
        end
        children = se_array(children_index);
        fprintf('\n## Super UC: %d\n',iFather);
        if sum(tree_model.brothers_mx(iFather,:)) > 0
            brother_id = find(tree_model.brothers_mx(iFather,:) == 1);
            [se_array{iFather},children] = se_array{iFather}.find_map(children, se_array{brother_id});
            se_array{iFather} = se_array{iFather}.assembly(children, se_array{brother_id});
            se_array{iFather} = se_array{iFather}.compute_modes(se_array{brother_id}, options.exppt);
            se_array{iFather} = se_array{iFather}.projection(se_array{brother_id});
        else
            [se_array{iFather},children] = se_array{iFather}.find_map(children, []);
            se_array{iFather} = se_array{iFather}.assembly(children, []);
            se_array{iFather} = se_array{iFather}.compute_modes([], options.exppt);
            se_array{iFather} = se_array{iFather}.projection([]);
        end
        se_array(children_index) = children;
    end
end

%% Output
frf = se_array{1}.get_frf(options);
