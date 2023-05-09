function [dofs_new, permuter] = calculate_dof_indices(Nodes,n_dof)
% Collect interior dofs in I
for i=1:length(Nodes.I)
    Dofs.I((i-1)*n_dof+1:i*n_dof)=(Nodes.I(i)*n_dof-n_dof+1):1:(Nodes.I(i)*n_dof);
end
% Separate boundary dofs
for i=1:length(Nodes.L)
    Dofs.L((i-1)*n_dof+1:i*n_dof)=(Nodes.L(i)*n_dof-n_dof+1):1:(Nodes.L(i)*n_dof);
end
for i=1:length(Nodes.R)
    Dofs.R((i-1)*n_dof+1:i*n_dof)=(Nodes.R(i)*n_dof-n_dof+1):1:(Nodes.R(i)*n_dof);
end
for i=1:length(Nodes.B)
    Dofs.B((i-1)*n_dof+1:i*n_dof)=(Nodes.B(i)*n_dof-n_dof+1):1:(Nodes.B(i)*n_dof);
end
for i=1:length(Nodes.T)
    Dofs.T((i-1)*n_dof+1:i*n_dof)=(Nodes.T(i)*n_dof-n_dof+1):1:(Nodes.T(i)*n_dof);
end
for i=1:length(Nodes.BL)
    Dofs.BL((i-1)*n_dof+1:i*n_dof)=(Nodes.BL(i)*n_dof-n_dof+1):1:(Nodes.BL(i)*n_dof);
end
for i=1:length(Nodes.BR)
    Dofs.BR((i-1)*n_dof+1:i*n_dof)=(Nodes.BR(i)*n_dof-n_dof+1):1:(Nodes.BR(i)*n_dof);
end
for i=1:length(Nodes.TL)
    Dofs.TL((i-1)*n_dof+1:i*n_dof)=(Nodes.TL(i)*n_dof-n_dof+1):1:(Nodes.TL(i)*n_dof);
end
for i=1:length(Nodes.TR)
    Dofs.TR((i-1)*n_dof+1:i*n_dof)=(Nodes.TR(i)*n_dof-n_dof+1):1:(Nodes.TR(i)*n_dof);
end
% Collect all boundary dofs in A
Dofs.A = [Dofs.BL,Dofs.B,Dofs.BR,Dofs.R,...
          Dofs.TR,Dofs.T,Dofs.TL,Dofs.L];
% Dofs.A = [Dofs.L,Dofs.R,Dofs.B,Dofs.T,Dofs.BL,Dofs.BR,Dofs.TR,Dofs.TL];
dofs_new.BL = [1:numel(Dofs.BL)]; 
dofs_new.B = dofs_new.BL(end) + [1:numel(Dofs.B)];
dofs_new.BR = dofs_new.B(end) + [1:numel(Dofs.BR)];
dofs_new.R = dofs_new.BR(end) + [1:numel(Dofs.R)];
dofs_new.TR = dofs_new.R(end) + [1:numel(Dofs.TR)];
dofs_new.T = dofs_new.TR(end) + [1:numel(Dofs.T)];
dofs_new.TL = dofs_new.T(end) + [1:numel(Dofs.TL)];
dofs_new.L = dofs_new.TL(end) + [1:numel(Dofs.L)];
dofs_new.I = dofs_new.L(end) + [1:numel(Dofs.I)];
dofs_new.A = [dofs_new.BL, dofs_new.B, dofs_new.BR, dofs_new.R,...
              dofs_new.TR, dofs_new.T, dofs_new.TL, dofs_new.L];
dofs_new.nDOF = numel([dofs_new.A,dofs_new.I]);
dofs_new.nDOF_A = numel([dofs_new.A]);

permuter([Dofs.A,Dofs.I],1:dofs_new.nDOF) = speye(dofs_new.nDOF);

end