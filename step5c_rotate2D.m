function [mid_profile_xy1,body_center_xy1,throat_center_xy1,surface_profile_xy1,body1_profile_xy1,body2_profile_xy1,angle_rotation]=step5c_rotate2D(sig1,mid_profile_xy,body_center_xy,throat_center_xy,surface_profile_xy,body1_profile_xy,body2_profile_xy,vis)
%This function rotates all 2D projected points such that body centers are horizontal. This is to ensure rotational invariance of the algorithm. It is part of Step 5c of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
X=body_center_xy(knnsearch(body_center_xy(:,2),max(body_center_xy(:,2))),1)-body_center_xy(knnsearch(body_center_xy(:,2),min(body_center_xy(:,2))),1);
Y=body_center_xy(knnsearch(body_center_xy(:,2),max(body_center_xy(:,2))),2)-body_center_xy(knnsearch(body_center_xy(:,2),min(body_center_xy(:,2))),2);
angle_rotation= atan2(norm(cross([X,Y,0],[1,0,0])),dot([X,Y,0],[1,0,0]));%Y is ensured to be positive, tan inverse is in the range o to pi in anticlockwise direction, correspondingly the rotation of the same angle is in clockwise direction
mid_profile_xy1=mid_profile_xy;body_center_xy1=body_center_xy;surface_profile_xy1=surface_profile_xy;body1_profile_xy1=body1_profile_xy;body2_profile_xy1=body2_profile_xy;

if angle_rotation~=0
    center = repmat([throat_center_xy(:,1); throat_center_xy(:,2)], 1, size(mid_profile_xy,1)).';
    theta = (round((angle_rotation)*10^(sig1)))/10^(sig1);% theta degrees rotation
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];vo = (mid_profile_xy - center)*R + center; % shift points so that center of rotation is at the origin,apply the rotation about the origin,shift origin back
    mid_profile_xy1(:,1) = vo(:,1);mid_profile_xy1(:,2) = vo(:,2);
    
    center = repmat([throat_center_xy(:,1); throat_center_xy(:,2)], 1, size(surface_profile_xy,1)).';
    vo = (surface_profile_xy - center)*R + center; % shift points so that center of rotation is at the origin,apply the rotation about the origin,shift origin back
    surface_profile_xy1(:,1) = vo(:,1);surface_profile_xy1(:,2) = vo(:,2);
    
    center = repmat([throat_center_xy(:,1); throat_center_xy(:,2)], 1, size(body_center_xy,1)).';
    vo = (body_center_xy - center)*R + center; % shift points so that center of rotation is at the origin,apply the rotation about the origin,shift origin back for body centers
    body_center_xy1(:,1) = vo(:,1);body_center_xy1(:,2) = vo(:,2);
    
    center = repmat([throat_center_xy(:,1); throat_center_xy(:,2)], 1, size(body1_profile_xy,1)).';
    vo = (body1_profile_xy - center)*R + center; % shift points so that center of rotation is at the origin,apply the rotation about the origin,shift origin back for body profile1
    body1_profile_xy1(:,1) = vo(:,1);body1_profile_xy1(:,2) = vo(:,2);
    
    center = repmat([throat_center_xy(:,1); throat_center_xy(:,2)], 1, size(body2_profile_xy,1)).';
    vo = (body2_profile_xy - center)*R + center; % shift points so that center of rotation is at the origin,apply the rotation about the origin,shift origin back for body profile 2
    body2_profile_xy1(:,1) = vo(:,1);body2_profile_xy1(:,2) = vo(:,2);
else
    mid_profile_xy1=mid_profile_xy;body_center_xy1=body_center_xy;surface_profile_xy1=surface_profile_xy;body1_profile_xy1=body1_profile_xy;body2_profile_xy1=body2_profile_xy;
end
throat_center_xy1= [mean(mid_profile_xy1(:,1)),mean(mid_profile_xy1(:,2))];

if vis ==1
    figure, scatter(mid_profile_xy1(:,1),mid_profile_xy1(:,2),'b','filled');hold on
    scatter(body_center_xy1(:,1),body_center_xy1(:,2),'r','filled');hold on
    scatter(throat_center_xy1(:,1),throat_center_xy1(:,2),'r','filled');hold on
    scatter(surface_profile_xy1(:,1),surface_profile_xy1(:,2),'y');hold on
    scatter(body1_profile_xy1(:,1),body1_profile_xy1(:,2),'g');hold on;
    scatter(body2_profile_xy1(:,1),body2_profile_xy1(:,2),'g');hold on;
    title('projected and rotated points on best fit 2D plane')
end
%---------------------------------------------------------------------------- END OF CODE ----------------------------------------------------------------------------------
