function [K_UC,M_UC,C_UC,UC_nodes,UC_coordinates,L_UC_x,L_UC_y] = import_FE_UC(filename,import_option)

switch import_option
    case 'mat'
        % Load the matrices
        load([filename,'.mat']);
        
        M_UC = model.M;
        K_UC = model.K+1i*model.K4gg;
        C_UC = model.C;
        
        % Symmetrize matrices (round-offs)
        M_UC = (M_UC+M_UC.')/2;
        K_UC = (K_UC+K_UC.')/2;
        C_UC = (C_UC+C_UC.')/2;
        
        % Load the coordinates
        UC_coordinates = model.coordinates;
    case 'op4'
        
        % Generate and load the matrices
        [~] = system(['"C:\Program Files\op4reader\op4reader.exe" modelmatrices.op4']);
        ModelMatrices = load('modelmatrices.mat');
        if isempty(ModelMatrices.Mgg)
            warning('!!! Empty mass matrix, check model.')
        elseif isempty(ModelMatrices.Kgg)
            warning('!!! Empty stiffness matrix, check model.')
        end
        
        % Symmetrize matrices
        M_UC = (ModelMatrices.Mgg+ModelMatrices.Mgg.')/2;
        K_UC = (ModelMatrices.Kgg+ModelMatrices.Kgg.')/2;
        C_UC = sparse(length(M_UC),length(M_UC));
        
        % Read model coordinates, node groups and DOF indices
        UC_coordinates = readcoords([filename,'.dat']);
        % Determine origin: unit cell is assumed to be oriented in the xy-plane
        direction.L1 = '+x'; % Left to Right
        direction.L2 = '+y'; % Top to Bottom
        direction.L3 = '+z';
        [~,iO]=intersect(UC_coordinates(:,2:3), min(UC_coordinates(:,2:3)),'rows');
        origin = UC_coordinates(iO,1);
        % Node groups
        Nodes = calculate_nodenumbers(direction,origin,UC_coordinates);
        
    otherwise
        error('Select right import_option')
end


UC_nodes.I = Nodes.interior;
UC_nodes.L = Nodes.left;
UC_nodes.R = Nodes.right;
UC_nodes.B = Nodes.bottom;
UC_nodes.T = Nodes.top;
UC_nodes.BL = Nodes.bottomleft;
UC_nodes.BR = Nodes.bottomright;
UC_nodes.TL = Nodes.topleft;
UC_nodes.TR = Nodes.topright;

L_UC_x = max(UC_coordinates(:,2));
L_UC_y = max(UC_coordinates(:,3));