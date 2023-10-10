function [neighbour]= step2_neighbours(separator,farthest_voxel,locals,voxels)
% This function finds the neighbours and coordination number of each pore volume to facilitate finding constriction surfaces between adjacent pore
% volumes. This is part of Step 2 of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%----------------------------------------------------------------- START CODE--------------------------------------------------------------------------
% The fartest pore volume voxel from the center of that pore volume is used to construct a sphere and all other pores inside that sphere are
% assigned as the neighbour of that pore

% figure,scatter3(voxels(locals,1),voxels(locals,2),voxels(locals,3),'filled','y');hold on;
% txt=num2str(linspace(1,size(locals,1),size(locals,1)).');
% text(voxels(locals,1),voxels(locals,2),voxels(locals,3),txt,'Color','red','FontSize',14);

neighbour=[];%coord_num=[];
%figure,
for i=1:size(locals,1)
    centerX = voxels(locals(i,1),1);centerY = voxels(locals(i,1),2);centerZ = voxels(locals(i,1),3);radius = farthest_voxel(i,1);
    sphereVoxels = (voxels(:,1) - centerX).^2 + (voxels(:,2) - centerY).^2 + (voxels(:,3) - centerZ).^2 <= radius.^2;
    sphereVoxels_index=unique(separator(sphereVoxels));% finding index within sphere
    k = arrayfun(@(a) (pdist2(voxels(locals(a),:),voxels(locals(i),:))),sphereVoxels_index(sphereVoxels_index>i));%distance between neighbouring pore centers
    neighbour=[neighbour; repmat(i,size(sphereVoxels_index(sphereVoxels_index>i),1),1) sphereVoxels_index(sphereVoxels_index>i)];% the current pore index, all neighbouring indices greater than the current index, distance between the neighbouring indices
    %coord_num=[coord_num; i, size(sphereVoxels_index,1)-1];
    % scatter3(voxels(sphereVoxels,1),voxels(sphereVoxels,2),voxels(sphereVoxels,3));hold on
    % pore_voxel= [voxels((separator==i),1) voxels((separator==i),2) voxels((separator==i),3)];
    % scatter3(pore_voxel(:,1),pore_voxel(:,2),pore_voxel(:,3),'r');hold on;
    %disp(['Portion of ',num2str(i),' out of ',num2str(length(locals)), 'processed']);
end
%----------------------------------------------------------------- END CODE--------------------------------------------------------------------------