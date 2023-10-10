function [points1,combined_profile_i]=step5i_inflexion3d(angle_rotation,throat_center_xy,points,mid_profile_xy1,mid_profile,body1_profile,body2_profile,mid_profile_xy,body1_profile_xy,body2_profile_xy,sig)
%This function finds 4 inflexion points backprojected in three dimensions. It is part of Step 5i of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
throat_axis_center=[mean(mid_profile_xy1(:,1)),mean(mid_profile_xy1(:,2))];
if angle_rotation ~= 0
    center = repmat([throat_center_xy(:,1); throat_center_xy(:,2)], 1, 4).';theta = (round((2*pi-angle_rotation)*10^(sig)))/10^(sig);% 270 degrees rotation
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];vo = (points(:,1:2) - center)*R + center; % shift points so that center of rotation is at the origin,apply the rotation about the origin,shift origin back
    points(:,1) = vo(:,1);points(:,2) = vo(:,2);
    
    center = throat_center_xy;theta = (round((2*pi-angle_rotation)*10^(sig)))/10^(sig);% 270 degrees rotation
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];vo = (throat_axis_center(:,1:2) - center)*R + center; % shift points so that center of rotation is at the origin,apply the rotation about the origin,shift origin back
    throat_axis_center(:,1) = vo(:,1);throat_axis_center(:,2) = vo(:,2);throat_axis_center=[throat_axis_center 0];
end
throat_axis_center=[throat_axis_center 0];

% 1st back projection to mid throat plane and then to voxels
combined_profile_xy= [mid_profile_xy;body1_profile_xy;body2_profile_xy];
combined_profile_i=[mid_profile;body1_profile;body2_profile];
for i=1:size(points,1)
    if points(i,3)==0
        IDX = knnsearch(mid_profile_xy,points(i,1:2));
        points1(i,:)= mid_profile(IDX,:);
    else
        IDX = knnsearch(combined_profile_xy,points(i,1:2));
        points1(i,:)= combined_profile_i(IDX,:);
    end
end

if throat_axis_center(1,3)==0
    IDX = knnsearch(mid_profile_xy,throat_axis_center(1,1:2));
    throat_axis_center= mid_profile(IDX,:);
else
    IDX = knnsearch(combined_profile_xy,throat_axis_center(1,1:2));
    throat_axis_center= combined_profile_i(IDX,:);
end
%--------------------------------------------------------------------------END CODE------------------------------------------------------------------------------------