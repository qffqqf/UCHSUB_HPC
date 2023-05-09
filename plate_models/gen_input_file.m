clc
clear all

nx = 30;
ny = 30;
pattern_uc = reshape([1:nx*ny], [nx,ny]);
pattern_bc = zeros([2*nx+1,2*ny+1]);
pattern_bc(:,1) = 1;
pattern_bc(:,end) = 1;
pattern_bc(1,:) = 1;
pattern_bc(end,:) = 1;
save('plate_30x30.mat','pattern_uc','pattern_bc')