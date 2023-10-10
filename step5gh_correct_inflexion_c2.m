function point= step5gh_correct_inflexion_c2(qci,tolc,name)
% This function finds point pertaining to condition 2 in correct_inflexion_points.m: find the leftmost point that slopes upward in 1st quadrant, for example
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
if (strcmp(name,'first')==1) || (strcmp(name,'fourth')==1)
    qci= sortrows(qci,1);
else
    qci= sortrows(qci,1,'descend');
end

if (strcmp(name,'first')==1) || (strcmp(name,'third')==1)
    for i=1:size(qci,1)
        if i== size(qci,1)
            point=[qci(i,1:2) 1];%scatter(point1(1, 1),point1(1, 2), 'm','filled');hold on;
            break;
        end
        if abs(qci(i+1,1)-qci(i,1))>=tolc && (qci(i+1,2)-qci(i,2))/(qci(i+1,1)-qci(i,1))>0
            point=[qci(i,1:2) 1];%scatter(point1(1, 1),point1(1, 2), 'm','filled');hold on;
            break
        else
            continue
        end
    end
end

if (strcmp(name,'second')==1) || (strcmp(name,'fourth')==1)
    for i=1:size(qci,1)
        if i== size(qci,1)
            point=[qci(i,1:2) 1];%scatter(point1(1, 1),point1(1, 2), 'm','filled');hold on;
            break;
        end
        if abs(qci(i+1,1)-qci(i,1))>=tolc && (qci(i+1,2)-qci(i,2))/(qci(i+1,1)-qci(i,1))<0
            point=[qci(i,1:2) 1];%scatter(point1(1, 1),point1(1, 2), 'm','filled');hold on;
            break
        else
            continue
        end
        if size(qci,1)==1
            point=[qci(1,1:2) 1];%scatter(point1(1, 1),point1(1, 2), 'm','filled');hold on;
        end
    end
end
%------------------------------------------------------------------- END CODE-------------------------------------------------------------------------------

