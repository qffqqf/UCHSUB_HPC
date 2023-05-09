clear all
clc
addpath('C:\Research\Projects\Substructuring\HSUB\Functions')
addpath('C:\Research\Projects\Substructuring\HSUB\plate_models\')
warning('off','all')
tic
%% UC model setup
filename = 'uclvb_spcoarse';
plate_file = 'plate_15x15.mat';
savename = "..\input_model_normal\input_model_15x15.mat";
import_option = 'mat'; 
n_dofnodes = 3; 
[K_UC_FOM, M_UC_FOM, C_UC_FOM, UC_nodes, UC_coordinates, L_UC_x, L_UC_y] = import_FE_UC3(filename);
load([filename,'.mat'])
%% Single UC model  
disp('Loading UC model ...')
[UC_dofs, permuter] = calculate_dof_indices(UC_nodes,n_dofnodes);

%% Get force
F_UC_FOM = zeros(UC_dofs.nDOF,1);
UC_forc = 6*(15-1)+1; 
picker_F = @(r) norm(r-[L_UC_x,L_UC_y,0])<1e-10;
force_UC_node = pick_nodes(mesh_data, picker_F);
F_UC_FOM(force_UC_node*3) = 1;

%% Get response
R_UC_FOM = zeros(UC_dofs.nDOF,1);
UC_resp = 10*(15-1)+1; 
picker_R = @(r) norm(r-[L_UC_x,L_UC_y,0])<1e-10;
response_UC_node = pick_nodes(mesh_data, picker_R);
R_UC_FOM(response_UC_node*3) = 1;

%% Permutation
uc_model.K_UC = permuter'*K_UC_FOM*permuter;
uc_model.M_UC = permuter'*M_UC_FOM*permuter;
uc_model.F_UC = permuter'*F_UC_FOM;
uc_model.R_UC = permuter'*R_UC_FOM;
uc_model.UC_dofs = UC_dofs;

%% Load plate model 
load(plate_file)
pattern_dofs = compute_connectivity_relations(pattern_uc);
assembly_info = get_assembly_info(UC_dofs, pattern_dofs, pattern_uc);
[uc_blocks, root_index, pointer_mx, brothers_mx] = get_quadtree(pattern_uc);
tree_model.assembly_info = assembly_info;
tree_model.uc_index = cell2mat(uc_blocks);
tree_model.root_index = root_index;
tree_model.pointer_mx = pointer_mx;
tree_model.brothers_mx = brothers_mx;

%% Parameters
param.forc = zeros(numel(pattern_uc),1);
param.forc(UC_forc) = 1;
param.resp = zeros(numel(pattern_uc),1);
param.resp(UC_resp) = 1;

%%  Save
save(savename, "uc_model", "tree_model", "param")
