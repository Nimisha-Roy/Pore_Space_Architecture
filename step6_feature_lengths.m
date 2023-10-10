function [size_body,size_throat]=step6_feature_lengths(throat_body_new2, num_pore_body2,voxels2,X_start1,X_end1, Y_start1,Y_end1, Z_start1,Z_end1,indices2)
%This function computes body and throat lengths in a,b,c axes of best ellipsoid fit. It excludes the features along the boundary of the structure. This is a postprocessing step of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
%----------------------------------------------------------------- START CODE--------------------------------------------------------------------------
throat_body_new2 = reshape(throat_body_new2, numel(throat_body_new2),1);
size_body=[];size_throat=[];
%pore_body_volume= zeros(length(throat_body_new),1); % to store throat_body values to compute individual pore volume later
%figure,
for u = 1:num_pore_body2 % to find volume of pore bodies by counting number of voxels for bodies and assigning the other throat voxels as zeros
    pos=(throat_body_new2==indices2(u));
    %check if pores are not along boundary and their sizes are atleast 10 voxels
    if sum(ismember(X_start1,voxels2(pos,1)))==0 && sum(ismember(X_end1,voxels2(pos,1)))==0 && sum(ismember(Y_start1,voxels2(pos,2)))==0 && sum(ismember(Y_end1,voxels2(pos,2)))==0 && sum(ismember(Z_start1,voxels2(pos,3)))==0 && sum(ismember(Z_end1,voxels2(pos,3)))==0 && sum(pos)>10
        ell = step6_ellipsoidfit(voxels2(throat_body_new2==indices2(u),:));
        %scatter3(voxels2(throat_body_new2==indices2(u),1),voxels2(throat_body_new2==indices2(u),2),voxels2(throat_body_new2==indices2(u),3),'MarkerFaceColor',rand(1,3));hold on
        size_body=[size_body;ell(1,4)*2, ell(1,5)*2,ell(1,6)*2];% all three axes length of ellipsoid
    end
end

%disp(['Portion of ',num2str(u),' out of ',num2str(num_pore_body), 'processed']);
k=1;
%figure,
for u = num_pore_body2+1:length(indices2)% to find volume of pore throats
    pos=(throat_body_new2==indices2(u));
    if sum(ismember(X_start1,voxels2(pos,1)))==0 && sum(ismember(X_end1,voxels2(pos,1)))==0 && sum(ismember(Y_start1,voxels2(pos,2)))==0 && sum(ismember(Y_end1,voxels2(pos,2)))==0 && sum(ismember(Z_start1,voxels2(pos,3)))==0 && sum(ismember(Z_end1,voxels2(pos,3)))==0 && sum(pos)>10
        ell = step6_ellipsoidfit(voxels2(throat_body_new2==indices2(u),:));
        %scatter3(voxels2(throat_body_new2==indices2(u),1),voxels2(throat_body_new2==indices2(u),2),voxels2(throat_body_new2==indices2(u),3),'MarkerFaceColor',rand(1,3));hold on
        size_throat=[size_throat;ell(1,4)*2, ell(1,5)*2,ell(1,6)*2];
    end
    k=k+1;
    %disp(['Portion of ',num2str(u),' out of ',num2str(length(indices)), 'processed']);
end
%----------------------------------------------------------------- END CODE--------------------------------------------------------------------------