clc
clear
close all

N=129;% grid size
Npts=N-2; % interior grid
h=1;
[X,Y,Z] = ndgrid(1:N,1:N,1:N);

T = double(niftiread('Brain01_1002.nii.gz'));
R = double(niftiread('Brain04_1009.nii.gz'));

T_info = niftiinfo('Brain01_1002.nii.gz');
R_info = niftiinfo('Brain04_1009.nii.gz');

niftiwrite(T, 'Brain01_1002_local.nii.gz');
niftiwrite(R, 'Brain04_1009_local.nii.gz');

[mx,my,mz]=size(R);
mlyr=2;
tic
[pos_xkFTR,pos_ykFTR,pos_zkFTR, TtoR, TtoR_jump]=Img_Reg3D_LagMul_faithful_v3_mtlres_example2(T,R,N,mlyr);
toc
figure(1)
gridplot3D_framed_axisR(pos_xkFTR, pos_ykFTR, pos_zkFTR,1,(N-1)/2+16,1),grid on
% gridplot3D_framed_axis(pos_xkFTR, pos_ykFTR, pos_zkFTR,1,(N-1)/2-16,3),grid on

[phi_jd_jump, ~, ~, ~]=compute_JD_and_Curl3D(pos_xkFTR,pos_ykFTR,pos_zkFTR,1);

% TtoR = imresize3(TtoR, [mx,my,mz]);
% TtoR_jump = imresize3(TtoR_jump, [mx,my,mz]);

niftiwrite(TtoR,'Brain01_1002_to_1009.nii.gz');
niftiwrite(TtoR_jump,'Brain01_1002_to_1009_jump.nii.gz');
niftiwrite(phi_jd_jump,'Brain01_1002_to_1009_jump_JD.nii.gz');


% j1 d1 Jm0.76231 Jmc-0.0074424 Jmf-0.33579 ts8.8197e-05 r0.29108 e47 t47 lyr128
% 
% r_local =
% 
%     0.7439
% 
% 
% r_final =
% 
%     0.2703
% 
% Elapsed time is 504.810984 seconds.




