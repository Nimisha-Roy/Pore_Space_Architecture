function [separator,farthest_voxel,locals] = step1_separate_porevolumes(locals,voxels,vox3d)
% This function discretizes the input pore space into pore volumes based on the computed local peaks in the LocalMaxima function. It is part of Step 1 of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
x_voxel = voxels(:,1);y_voxel = voxels(:,2);z_voxel = voxels(:,3); % voxel is linearized matrix of all voxels and vox3d is voxelized structure of voxels incorporating particles and pores
voxelInfo = reshape(vox3d, numel(vox3d),1);
N_voxel = length(x_voxel);separator = zeros(N_voxel,1);
npp = 20; % number of processing voxels more than which cannot be handled by memory
counter = 0:npp:size(voxelInfo,1);
if(counter(end) ~= size(voxelInfo,1))
    counter = [counter size(voxelInfo,1)];
end
peak_voxels = [x_voxel(locals),y_voxel(locals),z_voxel(locals)]; % peak points

for j = 1:length(counter)-1
    % The below code calculate the distance of local peak points from certain number of voxels and finds what local cluster a voxel belongs to by finding the peak which is at the minimum distance to that voxel.
    % Note that when the minimum distance is zero, in fact it is a solid point and we should assign zero for that value
    voxels1 = [x_voxel(counter(j)+1:counter(j+1)),y_voxel(counter(j)+1:counter(j+1)),z_voxel(counter(j)+1:counter(j+1))]; % a group of raw voxels
    BinaryVoxelPortion = voxelInfo(counter(j)+1:counter(j+1))'; % the same group of binarized voxels
    DIST = pdist2(peak_voxels,voxels1); % distance between peak voxels and the group of voxels
    [~,ind] = min(DIST,[],1); % ind stores the index of the assigned nearest peak to each voxel
    ind(BinaryVoxelPortion == 0) = 0; % to assign zero index value to solid voxels
    separator(counter(j)+1:counter(j+1)) = ind';
end

% to store the farthest voxel belonging to each peak voxel from each peak voxel
farthest_voxel=[];recordi=[];
for i=1:size(locals,1)
    s=voxels((separator==i),:);
    if size(s,1)==0
        recordi=[recordi;i];
        farthest_voxel=[farthest_voxel;0];
        continue
    end
    x=pdist2(voxels(locals(i),:),s);
    farthest_voxel=[farthest_voxel;max(x)];
    %disp(['Portion of ',num2str(i),' out of ',num2str(length(locals)), 'processed']);
end
if size(recordi,1)>0
    locals(recordi)=[];
    farthest_voxel(recordi)=[];
    separator=make_labels_sequential(separator);
end
%----------------------------------------------------------------- END CODE--------------------------------------------------------------------------
