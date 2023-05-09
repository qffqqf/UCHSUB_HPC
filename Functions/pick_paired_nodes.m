function nodes_list2_new = pick_paired_nodes(mesh_data, nodes_list1, nodes_list2, paired_info)

nNode = numel(nodes_list1);
nodes_list2_new = [];

for iNode = 1:nNode
    for jNode = 1:nNode
        ndID1 = nodes_list1(iNode);
        ndID2 = nodes_list2(jNode);
        coord1 = mesh_data.nd_data(ndID1,:);
        coord2 = mesh_data.nd_data(ndID2,:);
        if norm(coord1(paired_info) - coord2(paired_info))<1e-5%3e-4
            nodes_list2_new(end+1) = ndID2;
        end
    end
end

nodes_list2_new = nodes_list2_new.';