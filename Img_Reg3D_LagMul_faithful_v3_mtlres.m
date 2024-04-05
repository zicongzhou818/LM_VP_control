function [phi_final1, phi_final2, phi_final3, T_local, T_final] = Img_Reg3D_LagMul_faithful_v3_mtlres(T,R,N,mlyr)
T_original=T;
[mx,my,mz]=size(R);
ssd_origin = sum(sum(sum((T-R).^2)));
% T_final, phi_comp1re,phi_comp2re,phi_comp3re,T_final2, U1, U2, U3
%% Image registration 
h=1;
lll=0;
ts_init=5.0e-4;
for i = mlyr:-1:1
    lyr=(N-1)/(2^(i-1));
    sz=lyr+1;
    omega_size = (sz)^3;
%     ts = 1e-4;
    jd_bd = 0.2-lll*0.04;
    jd_bdm = 0.8;
    ts_tol = 8.0e-5;% tstep_ratio: stop when tstep less than tstep_ratio
    ts_up = 1.03;% tstep_up: magnification factor when ssd decreases
    ts_dn = 0.97;% tstep_down: reduction factor when doesn't ssd decrease
    imax = 1.0e2-lll*40;% ite_max: preset maximum number of iterations
    r_tol = 8e-2;% ratiotol: preset ratio to stop 
    r = 1;
    %Compute initial ssd on uniform grid
    
    % iteration flag
    better = 1;% Better: indicate whether the ssd decreases(Better=1) or not 
    ti = 0; % ii: total iteration steps
    ei = 0;% iii: iteration steps that ssd decreases
    % set up initial values for Phi(x,y)=(x+u1(x,y),y+u2(x,y))
    [X,Y,Z] = ndgrid(1:sz,1:sz,1:sz);% (i,j)--->(pos_x(i,j),pos_y(i,j)) on (1 to m,1 to n)
    phi_m1 = X;
    phi_m2 = Y;
    phi_m3 = Z;
    phi_comp1 = X;
    phi_comp2 = Y;
    phi_comp3 = Z;
    phi_final1 = phi_comp1;
    phi_final2 = phi_comp2;
    phi_final3 = phi_comp3;
    %% (ps: the final phi initially is the same as given phi here because the current phi is id)
    T_old = imresize3(T, [sz sz sz]);
    R_mtlres = imresize3(R, [sz sz sz]);
    ts = ts_init;
    T_old = interpn(X,Y,Z,T_old,phi_final1,phi_final2,phi_final3,'makima');
    T_local=T_old;
    ssd_old = sum(sum(sum((T_old-R_mtlres).^2)));
    ssd_initial = ssd_old;
    real_jsc=0;
    real_dice=0;
    Jmin_comp=1;
    while  (ts>ts_tol) && (ti<imax) && (r>r_tol) && (Jmin_comp>0)
        ti=ti+1; % record the number of iterations
        %% Step-2
        if better  % --(PS: Better holds for the first step in "while" loop)
            %% 3
            [Gx, Gy, Gz]=compute_gradient_vec(X, Y, Z, T_local, R_mtlres, phi_m1,phi_m2,phi_m3, phi_comp1, phi_comp2, phi_comp3,omega_size,h,lll);
        end
        %% 1
        phi_m1(2:sz-1, 2:sz-1, 2:sz-1) = X(2:sz-1, 2:sz-1, 2:sz-1) - ts*Gx(2:sz-1, 2:sz-1, 2:sz-1);
        phi_m2(2:sz-1, 2:sz-1, 2:sz-1) = Y(2:sz-1, 2:sz-1, 2:sz-1) - ts*Gy(2:sz-1, 2:sz-1, 2:sz-1);
        phi_m3(2:sz-1, 2:sz-1, 2:sz-1) = Z(2:sz-1, 2:sz-1, 2:sz-1) - ts*Gz(2:sz-1, 2:sz-1, 2:sz-1);
        [phi_m1,phi_m2,phi_m3] = gridunfolding(phi_m1,phi_m2,phi_m3,jd_bdm,lyr);
        %% Step-6  % phi(x,y)=X+U
        phi_comp1n = interpn(X, Y, Z, phi_comp1, phi_m1, phi_m2, phi_m3, 'makima'); 
        phi_comp2n = interpn(X, Y, Z, phi_comp2, phi_m1, phi_m2, phi_m3, 'makima'); 
        phi_comp3n = interpn(X, Y, Z, phi_comp3, phi_m1, phi_m2, phi_m3, 'makima'); 
        [phi_comp1n,phi_comp2n,phi_comp3n] = gridunfolding(phi_comp1n,phi_comp2n,phi_comp3n,jd_bd,lyr);
        %% Step-7
        T_temp = interpn(X,Y,Z, T_local, phi_m1, phi_m2, phi_m3,'makima'); 
        T_final = interpn(X,Y,Z,T_old, phi_comp1n, phi_comp2n, phi_comp3n,'makima');
        %% Step-8 Key part for gradient descent
        ssd=sum(sum(sum((T_temp-R_mtlres).^2))); % use sums to approximate integrals
        r_rel=ssd/ssd_initial;
        r=sum(sum(sum((T_final-R_mtlres).^2)))/ssd_initial; 
        [jsc, dice, ~, Jmin, ~, ~] = evaluations3(T_temp, R_mtlres, phi_m1,phi_m2,phi_m3);
        [~, ~, ~, Jmin_comp, ~, ~] = evaluations3(T_final, R_mtlres, phi_comp1n,phi_comp2n,phi_comp3n);

        if (ssd<ssd_old)
            ei=ei+1;
            T_local=T_temp;
            ts=ts*ts_up;
            ssd_old=ssd;
            better=1;
    
            real_jsc=jsc;
            real_dice=dice;
            
            phi_comp1=phi_comp1n; 
            phi_comp2=phi_comp2n; 
            phi_comp3=phi_comp3n; 
        else
            better=0;
            ts=ts*ts_dn;
        end 

        display(['j',num2str(real_jsc),' d',num2str(real_dice),' Jm',num2str(Jmin),' Jmc',num2str(Jmin_comp),' ts',num2str(ts),' r',num2str(r),' rr',num2str(r_rel),' e',num2str(ei),' t',num2str(ti),' lyr',num2str(lyr)]);
    end
    lll=lll+1;
    if lll~=mlyr
        phi_comp1 = interpn(phi_comp1, 'makima');
        phi_comp2 = interpn(phi_comp2, 'makima');
        phi_comp3 = interpn(phi_comp3, 'makima');
        phi_comp1 = (phi_comp1-1)/(sz-1);
        phi_comp2 = (phi_comp2-1)/(sz-1);
        phi_comp3 = (phi_comp3-1)/(sz-1);
        lyrn = (N-1)/(2^(i-2));

        phi_comp1 = phi_comp1*lyrn+1;
        phi_comp2 = phi_comp2*lyrn+1;
        phi_comp3 = phi_comp3*lyrn+1;

        phi_final1 = interpn(phi_final1, 'makima');
        phi_final2 = interpn(phi_final2, 'makima');
        phi_final3 = interpn(phi_final3, 'makima');

        phi_final1 = (phi_final1-1)/(sz-1);
        phi_final2 = (phi_final2-1)/(sz-1);
        phi_final3 = (phi_final3-1)/(sz-1);

        phi_final1 = phi_final1*lyrn+1;
        phi_final2 = phi_final2*lyrn+1;
        phi_final3 = phi_final3*lyrn+1;

        szn=lyrn+1;
        [X,Y,Z] = ndgrid(1:szn,1:szn,1:szn);
        phi_final1 = interpn(X, Y, Z, phi_final1, phi_comp1, phi_comp2, phi_comp3, 'makima'); 
        phi_final2 = interpn(X, Y, Z, phi_final2, phi_comp1, phi_comp2, phi_comp3, 'makima'); 
        phi_final3 = interpn(X, Y, Z, phi_final3, phi_comp1, phi_comp2, phi_comp3, 'makima'); 

        [phi_final1, phi_final2, phi_final3] = gridunfolding(phi_final1, phi_final2, phi_final3,jd_bd, lyr);
        T_fcheck = imresize3(T, [szn,szn,szn]);
        R_fcheck = imresize3(R, [szn szn szn]);
        [~, ~, ~, Jmin_final, ~, ~] = evaluations3(T_fcheck, R_fcheck, phi_final1, phi_final2, phi_final3);
    else
        [X,Y,Z] = ndgrid(1:sz,1:sz,1:sz);
        phi_final1 = interpn(X, Y, Z, phi_final1, phi_comp1, phi_comp2, phi_comp3, 'makima'); 
        phi_final2 = interpn(X, Y, Z, phi_final2, phi_comp1, phi_comp2, phi_comp3, 'makima'); 
        phi_final3 = interpn(X, Y, Z, phi_final3, phi_comp1, phi_comp2, phi_comp3, 'makima'); 
    end
    display(['j',num2str(real_jsc),' d',num2str(real_dice),' Jm',num2str(Jmin),' Jmc',num2str(Jmin_comp),' Jmf',num2str(Jmin_final),' ts',num2str(ts),' r',num2str(r),' e',num2str(ei),' t',num2str(ti),' lyr',num2str(lyr)]);

end

T = imresize3(T_original,[szn,szn,szn]);
T_final = interpn(X, Y, Z, T, phi_final1, phi_final2, phi_final3, 'makima');
T_local = imresize3(T_local,[mx,my,mz]);
T_final = imresize3(T_final,[mx,my,mz]);
r_local = sum(sum(sum((T_local-R).^2)))/ssd_origin
r_final = sum(sum(sum((T_final-R).^2)))/ssd_origin
end
