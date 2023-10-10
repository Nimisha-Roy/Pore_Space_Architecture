function [mid_profile_xy,body_center_xy,throat_center_xy,surface_profile_xy,body1_profile_xy,body2_profile_xy]=step5c_bestfit2D(locals,throat_center,unit_norm,x1,voxels,mid_profile,surface_profile,body1_profile,body2_profile,vis)
% This function projects throat, adjacent bodies and throat_body_surface voxels from best fit plane in 3d to xy, yz or xz plane to do computations in 2D with reference axes. This is Step 5c of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
mid_profile_xy = mid_profile; body_center_xy=[voxels(locals(throat_center(x1,5)),:);voxels(locals(throat_center(x1,6)),:)];
throat_center_xy=throat_center(x1,1:3);surface_profile_xy = surface_profile;
body1_profile_xy=body1_profile;body2_profile_xy=body2_profile;

% the direction towards which the unit normal is highly directed is the direction not counted for the plane
[~,ind_xy]=max(abs(unit_norm(:,1)));
if (ind_xy ==3);mid_profile_xy=[mid_profile_xy(:,1) mid_profile_xy(:,2)];
    body_center_xy=[body_center_xy(:,1) body_center_xy(:,2)];
    throat_center_xy=[throat_center_xy(:,1) throat_center_xy(:,2)];
    surface_profile_xy=[surface_profile_xy(:,1) surface_profile_xy(:,2)];
    body1_profile_xy=[body1_profile_xy(:,1) body1_profile_xy(:,2)];
    body2_profile_xy=[body2_profile_xy(:,1) body2_profile_xy(:,2)];
else if (ind_xy==2); mid_profile_xy=[mid_profile_xy(:,1) mid_profile_xy(:,3)];
        body_center_xy=[body_center_xy(:,1) body_center_xy(:,3)];
        throat_center_xy=[throat_center_xy(:,1) throat_center_xy(:,3)];
        surface_profile_xy=[surface_profile_xy(:,1) surface_profile_xy(:,3)];
        body1_profile_xy=[body1_profile_xy(:,1) body1_profile_xy(:,3)];
        body2_profile_xy=[body2_profile_xy(:,1) body2_profile_xy(:,3)];
    else
        mid_profile_xy=[mid_profile_xy(:,2) mid_profile_xy(:,3)];
        body_center_xy=[body_center_xy(:,2) body_center_xy(:,3)];
        throat_center_xy=[throat_center_xy(:,2) throat_center_xy(:,3)];
        surface_profile_xy=[surface_profile_xy(:,2) surface_profile_xy(:,3)];
        body1_profile_xy=[body1_profile_xy(:,2) body1_profile_xy(:,3)];
        body2_profile_xy=[body2_profile_xy(:,2) body2_profile_xy(:,3)];
    end
end

% to also project body centers on xy plane and  find the projected throat centers on xy plane
throat_center_xy= [mean(mid_profile_xy(:,1)),mean(mid_profile_xy(:,2))];
if vis ==1
    figure, scatter(mid_profile_xy(:,1),mid_profile_xy(:,2),'b');hold on;
    scatter(surface_profile_xy(:,1),surface_profile_xy(:,2),'y');hold on;
    scatter(body_center_xy(:,1),body_center_xy(:,2),'r', 'filled');hold on
    scatter(throat_center_xy(:,1),throat_center_xy(:,2),'r','filled');hold on
    scatter(body1_profile_xy(:,1),body1_profile_xy(:,2),'g');hold on;scatter(body2_profile_xy(:,1),body2_profile_xy(:,2),'g');hold on;
    line([body_center_xy(1,1),body_center_xy(2,1)],[body_center_xy(1,2),body_center_xy(2,2)]);
    title('projection of throat and voxels on best fit 2D plane')
end
%------------------------------------------------------------------------ END CODE --------------------------------------------------------------------------------