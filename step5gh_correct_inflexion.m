function [point1,point2,point3,point4]=step5gh_correct_inflexion(sizy,sizx,firstqci,secondqci,thirdqci,fourthqci,firstqtall,secondqtall,thirdqtall,fourthqtall,firstqti,secondqti,thirdqti,fourthqti,firstqcall,secondqcall,thirdqcall,fourthqcall)
% This function  finds the correct 4 inflexion points. It is part of Steps 5g and h of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
%if both throat and combined profile have inflexion points in the right direction and they are close, throat profile inflexion point is selected.
%If not close, throat profile inflexion point has boundary effects and so combined profile inflexion point is selected if it is the right direction. If there is no combined
%profile inflexion point but there is throat profile inflexion point, throat profile inflexion point is selected. If none of the profiles have inflexion point,the highest/lowest points of the corresponding/neighbouring quadrant is considered
%Right direction means, for example, for 1st quadrant in throat: leftmost point sloping upwards; for 1st quadrant in combined: leftmost point of sloping upwards

tol=sizx*0.02;%tolerance for closeness of points in throat boundary
tol_prot=(sizy/sizx)*0.5;%tolerance for protrusion
point1=[];point2=[];point3=[];point4=[];

%INFLEXION POINT CORRECTION
%condition 1- throat inflexion point considered %condition 2- combined inflexion point considered

%1st quadrant
if size(firstqti,1)~=0 && size(firstqci,1)~=0 % when both throat and combined profile have inflexion points
    % find protrusion ratio and distance between inflexion points of both profile
    [protrusion_ratio,dx]=step5gh_prot_rat_dx(firstqti,firstqci,firstqcall,'first');
end
if (size(firstqti,1)~=0 && size(firstqci,1)~=0 && (dx<tol || protrusion_ratio>tol_prot)) || (size(firstqci,1)==0 && size(firstqti,1)~=0) % if inflexion points of both profiles match or there are no inflexion points in combined profile, we consider throat profile inflexion point
    point1=step5gh_correct_inflexion_c1(firstqti,tol,'first');%condition 1 activated
else
    if size(firstqci,1)~=0
        point1=step5gh_correct_inflexion_c2(firstqci,tol,'first');
        %2nd sanity check that the point1 on combined profile indeed has a low protrusion ratio and is far away from the throat inflexion
        if size(firstqti,1)~=0
            [protrusion_ratio,dx]=step5gh_prot_rat_dx(firstqti,point1,firstqcall,'first');
            if (dx<tol || protrusion_ratio>tol_prot)
                point1=step5gh_correct_inflexion_c1(firstqti,tol,'first');
            end
        end
    else
        if size(firstqtall,1)>0%if both profiles don't have an inflexion point, we just consider the highest point of 1st quadrant
            [~,idx]= max(firstqtall(:,2));
            point1=[firstqtall(idx,:) 1];%scatter(point1(1, 1),point1(1, 2), 'm','filled');hold on;
        else
            if size(fourthqtall,1)>0
                [~,idx]= max(fourthqtall(:,2));
                point1=[fourthqtall(idx,:) 1];%scatter(point1(1, 1),point1(1, 2), 'm','filled');hold on;
            end
        end
    end
end

%2nd quadrant
if size(secondqti,1)~=0 && size(secondqci,1)~=0 % when both throat and combined profile have inflexion points
    % find protrusion ration and distance between inflexion points of both profile
    [protrusion_ratio,dx]=step5gh_prot_rat_dx(secondqti,secondqci,secondqcall,'second');
end
if (size(secondqti,1)~=0 && size(secondqci,1)~=0 && (dx<tol || protrusion_ratio>tol_prot)) || (size(secondqci,1)==0 && size(secondqti,1)~=0) % if inflexion points of both profiles match or there are no inflexion points in combined profile, we consider throat profile inflexion point
    point2=step5gh_correct_inflexion_c1(secondqti,tol,'second');%condition 1 activated
else if size(secondqci,1)~=0
        point2=step5gh_correct_inflexion_c2(secondqci,tol,'second');
        %2nd sanity check that the point1 on combined profile indeed has a low protrusion ratio and is far away from the throat inflexion
        if size(secondqti,1)~=0
            [protrusion_ratio,dx]=step5gh_prot_rat_dx(secondqti,point2,secondqcall,'second');
            if (dx<tol || protrusion_ratio>tol_prot)
                point2=step5gh_correct_inflexion_c1(secondqti,tol,'second');
            end
        end
    else if size(secondqtall,1)>0%if both profiles don't have an inflexion point, we just consider the highest point of 1st quadran
            [~,idx]= max(secondqtall(:,2));
            point2=[secondqtall(idx,:) 0];%scatter(point2(1, 1),point2(1, 2), 'm','filled');hold on;
        else if size(thirdqtall,1)>0
                [~,idx]= max(thirdqtall(:,2));
                point2=[thirdqtall(idx,:) 0];%scatter(point2(1, 1),point2(1, 2), 'm','filled');hold on;
            end
        end
        
    end
end

%3rd quadrant
if size(thirdqti,1)~=0 && size(thirdqci,1)~=0 % when both throat and combined profile have inflexion points
    % find protrusion ration and distance between inflexion points of both profile
    [protrusion_ratio,dx]=step5gh_prot_rat_dx(thirdqti,thirdqci,thirdqcall,'third');
end
if (size(thirdqti,1)~=0 && size(thirdqci,1)~=0 && (dx<tol || protrusion_ratio>tol_prot)) || (size(thirdqci,1)==0 && size(thirdqti,1)~=0) % if inflexion points of both profiles match or there are no inflexion points in combined profile, we consider throat profile inflexion point
    point3=step5gh_correct_inflexion_c1(thirdqti,tol,'third');%condition 1 activated
else if size(thirdqci,1)~=0
        point3=step5gh_correct_inflexion_c2(thirdqci,tol,'third');
        %2nd sanity check that the point1 on combined profile indeed has a low protrusion ratio and is far away from the throat inflexion
        if size(thirdqti,1)~=0
            [protrusion_ratio,dx]=step5gh_prot_rat_dx(thirdqti,point3,thirdqcall,'third');
            if (dx<tol || protrusion_ratio>tol_prot)
                point3=step5gh_correct_inflexion_c1(thirdqti,tol,'third');
            end
        end
    else if size(thirdqtall,1)>0
            [~,idx]= min(thirdqtall(:,2));
            point3=[thirdqtall(idx,:) 0];%scatter(point3(1, 1),point3(1, 2), 'm','filled');hold on;
        else if size(secondqtall,1)>0
                [~,idx]= min(secondqtall(:,2));
                point3=[secondqtall(idx,:) 0];%scatter(point3(1, 1),point3(1, 2), 'm','filled');hold on;
            end
        end
    end
end

%4th quadrant
if size(fourthqti,1)~=0 && size(fourthqci,1)~=0 % when both throat and combined profile have inflexion points
    % find protrusion ration and distance between inflexion points of both profile
    [protrusion_ratio,dx]=step5gh_prot_rat_dx(fourthqti,fourthqci,fourthqcall,'fourth');
end
if (size(fourthqti,1)~=0 && size(fourthqci,1)~=0 && (dx<tol || protrusion_ratio>tol_prot)) || (size(fourthqci,1)==0 && size(fourthqti,1)~=0) % if inflexion points of both profiles match or there are no inflexion points in combined profile, we consider throat profile inflexion point
    point4=step5gh_correct_inflexion_c1(fourthqti,tol,'fourth');%condition 1 activated
else if size(fourthqci,1)~=0
        point4=step5gh_correct_inflexion_c2(fourthqci,tol,'fourth');
        %2nd sanity check that the point1 on combined profile indeed has a low protrusion ratio and is far away from the throat inflexion
        if size(fourthqti,1)~=0
            [protrusion_ratio,dx]=step5gh_prot_rat_dx(fourthqti,point4,fourthqcall,'fourth');
            if (dx<tol || protrusion_ratio>tol_prot)
                point4=step5gh_correct_inflexion_c1(fourthqti,tol,'fourth');
            end
        end
    else  if size(fourthqtall,1)>0
            [~,idx]= min(fourthqtall(:,2));
            point4=[fourthqtall(idx,:) 0];%scatter(point4(1, 1),point4(1, 2), 'm','filled');hold on;
        else if size(firstqtall,1)>0
                [~,idx]= min(firstqtall(:,2));
                point4=[firstqtall(idx,:) 0];%scatter(point4(1, 1),point4(1, 2), 'm','filled');hold on;
            end
        end
    end
end

%--------------------------------------------------------------------------END CODE----------------------------------------------------------------------------------