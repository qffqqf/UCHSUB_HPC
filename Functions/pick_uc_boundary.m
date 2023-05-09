function [uc_nodes, uc_dofs] = pick_uc_boundary(mesh_data, uc_size)

eps = 1e-10;
Lx = uc_size(1);
Ly = uc_size(2);

%% pick functions
picker_i = @(r) (r(1)>0+eps) & (r(1)<Lx-eps) & (r(2)>0+eps) & (r(2)<Ly-eps);
picker_l = @(r) (abs(r(1)-0)<eps) & (abs(r(2)-0)*abs(r(2)-Ly)>eps);
picker_r = @(r) (abs(r(1)-Lx)<eps) & (abs(r(2)-0)*abs(r(2)-Ly)>eps);
picker_b = @(r) (abs(r(2)-0)<eps) & (abs(r(1)-0)*abs(r(1)-Lx)>eps);
picker_t = @(r) (abs(r(2)-Ly)<eps) & (abs(r(1)-0)*abs(r(1)-Lx)>eps);
picker_bl = @(r) norm([r(1)-0,r(2)-0])<eps;
picker_br = @(r) norm([r(1)-Lx,r(2)-0])<eps;
picker_tr = @(r) norm([r(1)-Lx,r(2)-Ly])<eps;
picker_tl = @(r) norm([r(1)-0,r(2)-Ly])<eps;

%% pick and pair
uc_nodes.I = pick_nodes(mesh_data, picker_i);
uc_nodes.L = pick_nodes(mesh_data, picker_l);
uc_nodes.R = pick_nodes(mesh_data, picker_r);
uc_nodes.B = pick_nodes(mesh_data, picker_b);
uc_nodes.T = pick_nodes(mesh_data, picker_t);
uc_nodes.BL = pick_nodes(mesh_data, picker_bl);
uc_nodes.BR = pick_nodes(mesh_data, picker_br);
uc_nodes.TR = pick_nodes(mesh_data, picker_tr);
uc_nodes.TL = pick_nodes(mesh_data, picker_tl);
uc_nodes.nNodes = size(mesh_data.nd_data, 1);

%% pair nodes
uc_nodes.R = pick_paired_nodes(mesh_data, uc_nodes.L, uc_nodes.R, [2,3]);
uc_nodes.T = pick_paired_nodes(mesh_data, uc_nodes.B, uc_nodes.T, [1,3]);
uc_nodes.BR = pick_paired_nodes(mesh_data, uc_nodes.BL, uc_nodes.BR, 3);
uc_nodes.TR = pick_paired_nodes(mesh_data, uc_nodes.BL, uc_nodes.TR, 3);
uc_nodes.TL = pick_paired_nodes(mesh_data, uc_nodes.BL, uc_nodes.TL, 3);

%% pick dofs
nd2dof = @(nd) nonzeros(reshape(mesh_data.nv_data(nd, :)',[],1));
uc_dofs.I = nd2dof(uc_nodes.I');
uc_dofs.L = nd2dof(uc_nodes.L');
uc_dofs.R = nd2dof(uc_nodes.R');
uc_dofs.B = nd2dof(uc_nodes.B');
uc_dofs.T = nd2dof(uc_nodes.T');
uc_dofs.BL = nd2dof(uc_nodes.BL');
uc_dofs.BR = nd2dof(uc_nodes.BR');
uc_dofs.TR = nd2dof(uc_nodes.TR');
uc_dofs.TL = nd2dof(uc_nodes.TL');
uc_dofs.A = [  uc_dofs.L; ...
               uc_dofs.R; ...
               uc_dofs.B; ...
               uc_dofs.T; ...
               uc_dofs.BL; ...
               uc_dofs.BR; ...
               uc_dofs.TR; ...
               uc_dofs.TL];
uc_dofs.nDOF = mesh_data.nDOF;         
uc_dofs.nBDOF = numel(uc_dofs.A);