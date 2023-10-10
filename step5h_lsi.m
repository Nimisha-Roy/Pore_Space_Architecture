function [ang1,ang2]=step5h_lsi(store_new,store_new1,points,vis)
% This function computes the longitudinal sharpness index of the throat boundary. It is part of Step 5h of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
point1=points(1,:);point2=points(2,:);point3=points(3,:);point4=points(4,:);
k= diff(store_new(:,2))./diff(store_new(:,1));
if point1(1,3)==0 && point2(1,3)==0
    A = [-k(knnsearch(store_new(:,1:2),point1(:,1:2))) 1; -k(knnsearch(store_new(:,1:2),point2(:,1:2))) 1];% linear equation of line at inflexion points tangential to the boundary at p1 and p2
    B = [ point1(1,2) - k(knnsearch(store_new(:,1:2),point1(:,1:2)))*point1(1,1) ; point2(1,2) - k(knnsearch(store_new(:,1:2),point2(:,1:2)))*point2(1,1)];
else
    if point1(1,3)==0 && point2(1,3)==1
        k1= diff(store_new1(:,2))./diff(store_new1(:,1));
        A = [-k(knnsearch(store_new(:,1:2),point1(:,1:2))) 1; -k1(knnsearch(store_new1(:,1:2),point2(:,1:2))) 1];% linear equation of line at inflexion points tangential to the boundary at p1 and p2
        B = [ point1(1,2) - k(knnsearch(store_new(:,1:2),point1(:,1:2)))*point1(1,1) ; point2(1,2) - k1(knnsearch(store_new1(:,1:2),point2(:,1:2)))*point2(1,1)];
    else
        if point1(1,3)==1 && point2(1,3)==0
            k1= diff(store_new1(:,2))./diff(store_new1(:,1));
            A = [-k1(knnsearch(store_new1(:,1:2),point1(:,1:2))) 1; -k(knnsearch(store_new(:,1:2),point2(:,1:2))) 1];% linear equation of line at inflexion points tangential to the boundary at p1 and p2
            B = [ point1(1,2) - k1(knnsearch(store_new1(:,1:2),point1(:,1:2)))*point1(1,1) ; point2(1,2) - k(knnsearch(store_new(:,1:2),point2(:,1:2)))*point2(1,1)];
        else
            k1= diff(store_new1(:,2))./diff(store_new1(:,1));
            A = [-k1(knnsearch(store_new1(:,1:2),point1(:,1:2))) 1; -k1(knnsearch(store_new1(:,1:2),point2(:,1:2))) 1];% linear equation of line at inflexion points tangential to the boundary at p1 and p2
            B = [ point1(1,2) - k1(knnsearch(store_new1(:,1:2),point1(:,1:2)))*point1(1,1) ; point2(1,2) - k1(knnsearch(store_new1(:,1:2),point2(:,1:2)))*point2(1,1)];
        end
    end
end
Xi = linsolve(A,B);x=[point1(1,1),Xi(1,1)];y=[point1(1,2), Xi(2,1)];
if vis ==1;plot(x,y, 'm', 'LineWidth',2);end
x=[point2(1,1),Xi(1,1)];y=[point2(1,2), Xi(2,1)];
if vis ==1;plot(x,y, 'm', 'LineWidth',2);end
DirVector1=[point1(1,1),point1(1,2)]-[Xi(1,1),Xi(2,1)];DirVector2=[point2(1,1),point2(1,1)]-[Xi(1,1),Xi(2,1)];
ang1=acosd( dot(DirVector1,DirVector2)/norm(DirVector1)/norm(DirVector2));
if point3(1,3)==0 && point4(1,3)==0
    A = [-k(knnsearch(store_new(:,1:2),point3(:,1:2))) 1; -k(knnsearch(store_new(:,1:2),point4(:,1:2))) 1];% linear equation of line at inflexion points tangential to the boundary at p3 and p4
    B = [ point3(1,2) - k(knnsearch(store_new(:,1:2),point3(:,1:2)))*point3(1,1) ; point4(1,2) - k(knnsearch(store_new(:,1:2),point4(:,1:2)))*point4(1,1)];
else
    if point3(1,3)==0 && point4(1,3)==1
        k1= diff(store_new1(:,2))./diff(store_new1(:,1));
        A = [-k(knnsearch(store_new(:,1:2),point3(:,1:2))) 1; -k1(knnsearch(store_new1(:,1:2),point4(:,1:2))) 1];% linear equation of line at inflexion points tangential to the boundary at p3 and p4
        B = [ point3(1,2) - k(knnsearch(store_new(:,1:2),point3(:,1:2)))*point3(1,1) ; point4(1,2) - k1(knnsearch(store_new1(:,1:2),point4(:,1:2)))*point4(1,1)];
    else
        if point3(1,3)==1 && point4(1,3)==0
            k1= diff(store_new1(:,2))./diff(store_new1(:,1));
            A = [-k1(knnsearch(store_new1(:,1:2),point3(:,1:2))) 1; -k(knnsearch(store_new(:,1:2),point4(:,1:2))) 1];% linear equation of line at inflexion points tangential to the boundary at p3 and p4
            B = [ point3(1,2) - k1(knnsearch(store_new1(:,1:2),point3(:,1:2)))*point3(1,1) ; point4(1,2) - k(knnsearch(store_new(:,1:2),point4(:,1:2)))*point4(1,1)];
        else
            k1= diff(store_new1(:,2))./diff(store_new1(:,1));
            A = [-k1(knnsearch(store_new1(:,1:2),point3(:,1:2))) 1; -k1(knnsearch(store_new1(:,1:2),point4(:,1:2))) 1];% linear equation of line at inflexion points tangential to the boundary at p3 and p4
            B = [ point3(1,2) - k1(knnsearch(store_new1(:,1:2),point3(:,1:2)))*point3(1,1) ; point4(1,2) - k1(knnsearch(store_new1(:,1:2),point4(:,1:2)))*point4(1,1)];
        end
    end
end

Xi = linsolve(A,B);x=[point3(1,1),Xi(1,1)];y=[point3(1,2), Xi(2,1)];
if vis==1;plot(x,y, 'm', 'LineWidth',2);end
x=[point4(1,1),Xi(1,1)];y=[point4(1,2), Xi(2,1)];
if vis==1;plot(x,y, 'm', 'LineWidth',2);end
DirVector1=[point3(1,1),point3(1,2)]-[Xi(1,1),Xi(2,1)];
DirVector2=[point4(1,1),point4(1,1)]-[Xi(1,1),Xi(2,1)];
ang2=acosd( dot(DirVector1,DirVector2)/norm(DirVector1)/norm(DirVector2));

%------------------------------------------------------------------------------END CODE-------------------------------------------------------------------