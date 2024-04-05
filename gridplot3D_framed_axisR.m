function gridplot3D_framed_axisR(X,Y,Z,h,k,orient)
% h2=2*h;
h2=h;
if orient==3
    [m,n,~] = size(X);
    for i = 1:h:m
        for j=1:h:n
            if i==1 || j==1 || mod(i,2)==0 || mod(j,h)
                a=squeeze(X(i, 1:h2:end, k));
                b=squeeze(Y(i, 1:h2:end, k));
                c=squeeze(Z(i, 1:h2:end, k));
                plot3(a, b, c,'r'), hold on;
                e=squeeze(X(1:h2:end, j, k));
                d=squeeze(Y(1:h2:end, j, k));
                f=squeeze(Z(1:h2:end, j, k));
                plot3(e, d, f,'r'), hold on;
            end
        end
    end
    axis([0, n, 0, n, 0, n])
elseif orient==2
    [m,~,n] = size(X);
    for i = 1:h:m
        for j=1:h:n
            if i==1 || j==1 || mod(i,2)==0 || mod(j,h)
                a=squeeze(X(i, k, 1:h2:end));
                b=squeeze(Y(i, k, 1:h2:end));
                c=squeeze(Z(i, k, 1:h2:end));
                plot3(a, b, c,'r'), hold on;
                e=squeeze(X(1:h2:end, k, j));
                d=squeeze(Y(1:h2:end, k, j));
                f=squeeze(Z(1:h2:end, k, j));
                plot3(e, d, f,'r'), hold on;
            end
        end
    end
    axis([0, n, 0, n, 0, n])    
elseif orient==1    
    [~,m,n] = size(X);
    for i = 1:h:m
        for j=1:h:n
            if i==1 || j==1 || mod(i,2)==0 || mod(j,h)
                a=squeeze(X(k, i, 1:h2:end));
                b=squeeze(Y(k, i, 1:h2:end));
                c=squeeze(Z(k, i, 1:h2:end));
                plot3(a, b, c,'r'), hold on;
                e=squeeze(X(k, 1:h2:end, j));
                d=squeeze(Y(k, 1:h2:end, j));
                f=squeeze(Z(k, 1:h2:end, j));
                plot3(e, d, f,'r'), hold on;
            end
        end
    end
    axis([0, n, 0, n, 0, n])      
end