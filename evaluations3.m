function [j, d, fmax, fmin, cv_max, cv_min] = evaluations3(T, R, Phi1,Phi2,Phi3)
    [f0, cv1,cv2,cv3] = compute_JD_and_Curl3D(Phi1,Phi2,Phi3,1);
    fmax=max(max(max(f0)));
    fmin=min(min(min(f0)));

    norm_cv=cv1.^2+cv2.^2+cv3.^2;
    cv_max=max(norm_cv,[],'all','omitnan');
    if cv_max>0
        norm_cv(norm_cv<=0)=nan;
    end
    cv_min=min(norm_cv,[],"all",'omitnan');

    thresh = graythresh(T);
    BW_T = imbinarize(T,thresh);
    thresh = graythresh(R);
    BW_R = imbinarize(R,thresh);

    j = jaccard(BW_T,BW_R);
    d = dice(BW_T,BW_R);
end

