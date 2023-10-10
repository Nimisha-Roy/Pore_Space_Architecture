function [protrusion_ratio,dx]=step5gh_prot_rat_dx(qti,qci,qcall,name)
%This function finds the protusion ratio and distance between line joining 2 points and boundary. This is part of Step 5g,h of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
if (strcmp(name,'first')==1) || (strcmp(name,'fourth')==1)
    qti= sortrows(qti,1);qci= sortrows(qci,1);
else
    qti= sortrows(qti,1,'descend');qci= sortrows(qci,1,'descend');
end
dx= abs(qti(1,1)-qci(1,1));%x distance between the leftmost firstq and firstq2
qcall=qcall((qcall(:,1)>=min(qti(1,1),qci(1,1))) & (qcall(:,1)<=max(qti(1,1),qci(1,1))),:);
a=[(qti(1,1:2)-qci(1,1:2)),0];
b=[unique(qcall-qci(1,1:2),'rows')];b=[b,zeros(size(b,1),1)];
yheights=arrayfun(@(rowid) norm(cross(a,b(rowid,:)))/norm(a),(1:size(b,1)).');%projection of smallest y point in 1st quadrant and inflexion points of throat profile and combined profile
protrusion_ratio=max(yheights)/pdist2(qti(1,1:2),qci(1,1:2));
%------------------------------------------------------------------------------END CODE------------------------------------------------------------------------------
