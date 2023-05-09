function [father_blocks, father_index, pointer_mx, brothers_mx] = get_quadtree(pattern_uc)

pointer.I = [];
pointer.J = [];
pointer.V = [];
brothers.I = [];
brothers.J = [];
brothers.V = [];
pointer.nSE = 1;
father_blocks = {pattern_uc};
father_index = [1];
all_blocks = {pattern_uc};
all_index = [1];
for iStep = 1:100
    [child_blocks, child_index, pointer] = recursive_decomp(father_blocks, father_index, pointer);
    all_blocks = [all_blocks, child_blocks];
    all_index = [all_index, child_index];
    if numel(father_blocks) ~= numel(child_blocks)
        father_blocks = child_blocks;
        father_index = child_index;
    else
        break;
    end
end
[all_index, non_rep] = unique(all_index);
all_blocks = all_blocks(non_rep);


black_list = [];
for iChild = numel(all_blocks):-1:1
    if sum(black_list==iChild) == 0
        for jChild = iChild-1:-1:1
            if norm(size(all_blocks{iChild}) - size(all_blocks{jChild})) < eps
                black_list(end+1) = jChild;
                brothers.I(end+1) = all_index(jChild);
                brothers.J(end+1) = all_index(iChild);
                brothers.V(end+1) = 1;
            end
        end
    end
end
pointer_mx = sparse(pointer.I, pointer.J, pointer.V, max(father_index), max(pointer.J));
brothers_mx = (sparse(brothers.I, brothers.J, brothers.V, max(father_index), max(brothers.J)) > 0);

repeat_ = sum(brothers_mx,2);
black_list = [];
for iRep = 1:numel(repeat_)
    if repeat_(iRep) > 1
        rep_ = find(brothers_mx(iRep, :)>0);
        for iBro = 1:numel(rep_)-1
            if sum(black_list==rep_(iBro)) == 0
                black_list(end+1) = rep_(iBro);
                brothers_mx(:,rep_(end)) = brothers_mx(:,rep_(end)) + brothers_mx(:,rep_(iBro));
                brothers_mx(:,rep_(iBro)) = 0;
            end
        end
    end
end


function blocks = decomp(block0)
nx = size(block0,1);
ny = size(block0,2);
if and(nx > 1, ny > 1)
    nx1 = [1:round(nx/2)];
    nx2 = [1+round(nx/2):nx];
    ny1 = [1:round(ny/2)];
    ny2 = [1+round(ny/2):ny];
    block1 = block0(nx2,ny1);
    block2 = block0(nx2,ny2);
    block3 = block0(nx1,ny2);
    block4 = block0(nx1,ny1);
    blocks = {block1,block2,block3,block4};
elseif and(nx > 1, ny == 1)
    nx1 = [1:round(nx/2)];
    nx2 = [1+round(nx/2):nx];
    block1 = block0(nx2,1);
    block2 = block0(nx1,1);
    blocks = {block1,block2};
elseif and(nx == 1, ny > 1)
    ny1 = [1:round(ny/2)];
    ny2 = [1+round(ny/2):ny];
    block1 = block0(1,ny1);
    block2 = block0(1,ny2);
    blocks = {block1,block2};
else
    disp("what the hell")
end
    
function flag = inspection(block0)
nx = size(block0,1);
ny = size(block0,2);
if or(nx > 1, ny > 1)
    flag = 1;
else
    flag = 0;
end


function [child_blocks, child_index, pointer] = recursive_decomp(father_blocks, father_index, pointer)
child_blocks = {};
child_index = [];
for iFather = 1:numel(father_blocks)
    if inspection(father_blocks{iFather})
        children = decomp(father_blocks{iFather});
        for iChild = 1:numel(children)
            child_blocks{end+1} = children{iChild};
            child_index(end+1) = pointer.nSE + 1;
            pointer.nSE = pointer.nSE + 1;
            pointer.I(end+1) = child_index(end);
            pointer.J(end+1) = father_index(iFather);
            pointer.V(end+1) = 1;
        end
    else
        child_blocks{end+1} = father_blocks{iFather};
        child_index(end+1) = father_index(iFather);
    end
end






