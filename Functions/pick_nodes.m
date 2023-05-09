function nodes_array = pick_nodes(mesh_data, picker_func)

nNode = size(mesh_data.nd_data, 1);

nodes_array = [];
for iNode = 1:nNode
    coord = mesh_data.nd_data(iNode,:);
    if picker_func(coord)
        nodes_array(end+1) = iNode;
    end
end

nodes_array = nodes_array.';