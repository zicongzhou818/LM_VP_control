function [Gx, Gy, Gz]=compute_gradient_vec(X_init, Y_init, Z_init, T_local, R_mtlres, phi_m1,phi_m2,phi_m3, phi_comp1, phi_comp2, phi_comp3,omega_size,h,lll)

T_temp = interpn(X_init,Y_init,Z_init,T_local,phi_m1,phi_m2,phi_m3,'makima'); 
TR_diff = T_temp-R_mtlres; 
[T_y, T_x, T_z] = gradient(T_temp); 

Ix = (TR_diff.*T_x)*0.1;
Iy = (TR_diff.*T_y)*0.1;
Iz = (TR_diff.*T_z)*0.1;
    
[phi_m1_y,phi_m1_x,phi_m1_z]=gradient(phi_m1,h);
[phi_m2_y,phi_m2_x,phi_m2_z]=gradient(phi_m2,h);
[phi_m3_y,phi_m3_x,phi_m3_z]=gradient(phi_m3,h);

[phi_comp_jd, phi_comp_cv1, phi_comp_cv2, phi_comp_cv3]=compute_JD_and_Curl3D(phi_comp1, phi_comp2, phi_comp3,h);
lambda1 = phi_comp_jd/(sum(sum(sum(phi_comp_jd))))*omega_size;
lambda2x = phi_comp_cv1;
lambda2y = phi_comp_cv2;
lambda2z = phi_comp_cv3;

iix_iiix1 = lambda1.*(phi_m2_y.*phi_m3_z - phi_m2_z.*phi_m3_y);
iix_iiix2 = lambda1.*(phi_m3_x.*phi_m2_z - phi_m3_z.*phi_m2_x) - lambda2z;
iix_iiix3 = lambda1.*(phi_m2_x.*phi_m3_y - phi_m2_y.*phi_m3_x) + lambda2y;
[~,iix_iiix1_x,~] = gradient(iix_iiix1,h);
[iix_iiix2_y,~,~] = gradient(iix_iiix2,h);
[~,~,iix_iiix3_z] = gradient(iix_iiix3,h);
IIx_IIIx = iix_iiix1_x + iix_iiix2_y + iix_iiix3_z;

iiy_iiiy1 = lambda1.*(phi_m1_z.*phi_m3_y - phi_m3_z.*phi_m2_y) + lambda2z;
iiy_iiiy2 = lambda1.*(phi_m1_x.*phi_m3_z - phi_m1_z.*phi_m2_z);
iiy_iiiy3 = lambda1.*(phi_m1_y.*phi_m3_x - phi_m1_x.*phi_m3_y) - lambda2x;
[~,iiy_iiiy1_x,~] = gradient(iiy_iiiy1,h);
[iiy_iiiy2_y,~,~] = gradient(iiy_iiiy2,h);
[~,~,iiy_iiiy3_z] = gradient(iiy_iiiy3,h);
IIy_IIIy = iiy_iiiy1_x + iiy_iiiy2_y + iiy_iiiy3_z;

iiz_iiiz1 = lambda1.*(phi_m1_y.*phi_m2_z - phi_m1_z.*phi_m2_y) - lambda2y;
iiz_iiiz2 = lambda1.*(phi_m1_z.*phi_m2_x - phi_m1_x.*phi_m2_z) + lambda2x;
iiz_iiiz3 = lambda1.*(phi_m1_x.*phi_m2_y - phi_m1_y.*phi_m2_x);
[~,iiz_iiiz1_x,~] = gradient(iiz_iiiz1,h);
[iiz_iiiz2_y,~,~] = gradient(iiz_iiiz2,h);
[~,~,iiz_iiiz3_z] = gradient(iiz_iiiz3,h);
IIz_IIIz = iiz_iiiz1_x + iiz_iiiz2_y + iiz_iiiz3_z;

% Gx = imgaussfilt3(Ix+ IIx_IIIx, exp(1)) ;
% Gy = imgaussfilt3(Iy+ IIy_IIIy, exp(1)) ;
% Gz = imgaussfilt3(Iz+ IIz_IIIz, exp(1)) ;

Gx = imgaussfilt3(Ix + IIx_IIIx, exp(1)-(lll*0.5));
Gy = imgaussfilt3(Iy + IIy_IIIy, exp(1)-(lll*0.5));
Gz = imgaussfilt3(Iz + IIz_IIIz, exp(1)-(lll*0.5));

end