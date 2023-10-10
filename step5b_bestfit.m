function [mid_profile,surface_profile,unit_norm,body1_profile,body2_profile]=step5b_bestfit(set_body1,set_body2,set2,ste,x1,set,voxels,locations,locals,throat_center,vis)
%This function finds the best fit plane approximating throat and body voxels and passing through the two body centers in the direction of the length.
% It returns the unit normal of the plane and the projected throat,body and surface profiles. It is part of Step 5b of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------ BEGIN CODE ---------------------------------------------------------------------------
% This is the Solution of constrained least squares problem solved using lagrangian function with 2 lagrangian multipliers
A=voxels(locations,:);b=ones(size(locations,1),1);C=[voxels(locals(throat_center(x1,5)),:);voxels(locals(throat_center(x1,6)),:)];d=ones(2,1);
frfc= 2*A.'*A; frsc=C.'; srfc= C; srsc=zeros(2); frfc1= 2*A.'*b; srfc1= d;
P=[frfc frsc; srfc srsc]; P_inv= inv(P);frfc= P_inv(1:3,1:3); frsc=P_inv(1:3,4:5); srfc= P_inv(4:5,1:3); srsc=P_inv(4:5,4:5);
coeff= frfc*frfc1 + frsc*srfc1; unit_norm= coeff/(coeff(1)^2+coeff(2)^2+coeff(3)^2); % coeff is the normal to the plane with d=-1
set1=[voxels(locations,:);voxels(locals(throat_center(x1,5)),:);voxels(locals(throat_center(x1,6)),:)];
[xx,yy] = ndgrid(-ceil(max(set1(:)))-100*ste:ceil(max(set1(:)))+100*ste); % Providing a gridwork spanning 10 voxels on either side of body centers
z= (-coeff(1)*xx - coeff(2)*yy + 1)/coeff(3);%surf(xx,yy,z);

% projection of throat points on the plane
pt= voxels(locals(throat_center(x1,5)),:); dot_P= dot((voxels(locations,:)-pt),repmat((unit_norm.'),size(locations,1),1),2);%(P-pt).n

% projected points of the 3D points obtained by projecting the vector joining 3D points and origin onto normal of the plane and subtracting it from the 3D points
mid_profile=((voxels(locations,:)-pt)-bsxfun(@times, dot_P, (coeff.')))+repmat(pt,size(locations,1),1);

% to find the projected throat-pore surface voxels
pt= voxels(locals(throat_center(x1,5)),:); dot_P= dot((set2-pt),repmat((unit_norm.'),size(set2,1),1),2);%(P-pt).n
surface_profile=((set2-pt)-bsxfun(@times, dot_P, (coeff.')))+repmat(pt,size(set2,1),1);

% to find projected voxels of the adjacent pore bodies
pt= voxels(locals(throat_center(x1,5)),:); dot_P= dot((set_body1-pt),repmat((unit_norm.'),size(set_body1,1),1),2);%(P-pt).n
body1_profile=((set_body1-pt)-bsxfun(@times, dot_P, (coeff.')))+repmat(pt,size(set_body1,1),1);
pt= voxels(locals(throat_center(x1,5)),:); dot_P= dot((set_body2-pt),repmat((unit_norm.'),size(set_body2,1),1),2);%(P-pt).n
body2_profile=((set_body2-pt)-bsxfun(@times, dot_P, (coeff.')))+repmat(pt,size(set_body2,1),1);

if vis ==1
    figure, scatter3(set(:,1),set(:,2), set(:,3)); hold on;scatter3(set2(:,1),set2(:,2), set2(:,3)); hold on
    scatter3(mid_profile(:,1),mid_profile(:,2),mid_profile(:,3),'g','filled');hold on
    scatter3(voxels(locals(throat_center(x1,5)),1),voxels(locals(throat_center(x1,5)),2), voxels(locals(throat_center(x1,5)),3),'r', 'filled');hold on
    scatter3(voxels(locals(throat_center(x1,6)),1),voxels(locals(throat_center(x1,6)),2), voxels(locals(throat_center(x1,6)),3), 'r', 'filled');hold on
    scatter3(throat_center(x1,1),throat_center(x1,2),throat_center(x1,3), 'r', 'filled');hold on
    title('projection of throat and body voxels on best fit 3D plane')
end
%------------------------------------------------------------- END CODE ------------------------------------------------------------------------------------