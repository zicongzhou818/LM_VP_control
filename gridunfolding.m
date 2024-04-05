function [phi_m1, phi_m2, phi_m3] = gridunfolding(phi_m1, phi_m2, phi_m3, j_bd, lyr)
[JDm, ~,~,~] = compute_JD_and_Curl3D(phi_m1,phi_m2,phi_m3,1);
linear_indices = find(JDm<j_bd);
if isempty(linear_indices)==0 
    disp('unfolding');
    [xi,yj,zk]=ind2sub(size(JDm), linear_indices);
    for ith = 1: length(xi)
        ii=xi(ith);
        jj=yj(ith);
        kk=zk(ith);
        if ii<lyr+1 && jj<lyr+1 && kk<lyr+1
            phi_m1(ii,jj,kk) = (phi_m1(ii-1,jj,kk) + phi_m1(ii+1,jj,kk))*0.5;
            phi_m2(ii,jj,kk) = (phi_m2(ii,jj-1,kk) + phi_m2(ii,jj+1,kk))*0.5;
            phi_m3(ii,jj,kk) = (phi_m3(ii,jj,kk-1) + phi_m3(ii,jj,kk+1))*0.5;
        end
    end
end
end
