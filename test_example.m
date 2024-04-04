clc
clear
close all

N=129;% grid size
Npts=N-2; % interior grid
h=1;
[X,Y,Z] = ndgrid(1:N,1:N,1:N);

T = niftiread('reg_T1_brain7.nii');
R = niftiread('reg_T1_brain14.nii');

T_info = niftiinfo('reg_T1_brain7.nii');
R_info = niftiinfo('reg_T1_brain14.nii');

niftiwrite(T,'test00_fullview_b7.nii',T_info)
niftiwrite(R,'test00_fullview_b14.nii',R_info)



mlyr=2;
tic
[pos_xkFTR,pos_ykFTR,pos_zkFTR, TtoR, TtoR_jump]=Img_Reg3D_LagMul_faithful_v3_mtlres(T,R,N,mlyr);
toc
figure(1)
gridplot3D_framed_axisR(pos_xkFTR, pos_ykFTR, pos_zkFTR,1,(N-1)/2+16,1),grid on
% gridplot3D_framed_axis(pos_xkFTR, pos_ykFTR, pos_zkFTR,1,(N-1)/2-16,3),grid on

[phi_jd_jump, ~, ~, ~]=compute_JD_and_Curl3D(pos_xkFTR,pos_ykFTR,pos_zkFTR,1);
niftiwrite(TtoR,'test01_fullview_b7Tob14.nii', R_info)
niftiwrite(TtoR_jump,'test02_fullview_b7Tob14_jump.nii', R_info)
niftiwrite(phi_jd_jump,'test03_fullview_b7Tob14_jump_JD.nii')


% j0.86664 d0.92855 Jm0.90672 Jmc-0.0055505 Jmf0.010343 ts0.0020059 r0.31743 e47 t47 lyr128
% 
% r_local =
% 
%     0.3287
% 
% 
% r_final =
% 
%     0.2892
% 
% Elapsed time is 535.773667 seconds.




