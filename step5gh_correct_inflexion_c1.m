function point=step5gh_correct_inflexion_c1(qti,tol,name)
%This function finds point pertaining to condition 1 in correct_inflexion_points.m:  find the rightmost point among the leftmost set of inflexion points in quadrant 1, for example
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
if (strcmp(name,'first')==1) || (strcmp(name,'fourth')==1)
    qti= sortrows(qti,1);
else
    qti= sortrows(qti,1,'descend');
end
for i=1:size(qti,1)
    if i== size(qti,1)
        point=[qti(i,:) 0];
        break;
    end
    if abs(qti(i+1,1)-qti(i,1))<tol
        continue;
    else
        point=[qti(i,:) 0];
        break;
    end
end
if size(qti,1)==1
    point=[qti(1,:) 0];
end
%------------------------------------------------------------------- END CODE-------------------------------------------------------------------------------
