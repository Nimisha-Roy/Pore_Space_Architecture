function [throat_body,throat_center,neighbour,locals]= step4_throat_body(locals,voxels,throat_center,vox3d,neighbour,num_pore_body)
%This function discretizes the pore structure into pore bodies and pore throats based on MEDA, updates body and throat centers along with their neighbours. It is step 4 of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
x_voxel = voxels(:,1);y_voxel = voxels(:,2);z_voxel = voxels(:,3); % voxel is linearized matrix of all voxels and vox3d is voxelized structure of voxels incorporating particles and pores
voxelInfo = reshape(vox3d, numel(vox3d),1);
N_voxel = length(x_voxel);throat_body = zeros(N_voxel,1);
npp = 20; % number of processing voxels more than which cannot be handled by memory
counter = 0:npp:size(voxelInfo,1);
if(counter(end) ~= size(voxelInfo,1))
    counter = [counter size(voxelInfo,1)];
end
peak_voxels = [x_voxel(locals),y_voxel(locals),z_voxel(locals);x_voxel(throat_center(:,4)),y_voxel(throat_center(:,4)),z_voxel(throat_center(:,4)) ]; % peak points

for j = 1:length(counter)-1
    % The below code calculate the distance of local peak points from certain number of voxels and finds what local cluster a voxel belongs to by finding the peak which is at the minimum distance to that voxel.
    % Note that when the minimum distance is zero, in fact it is a solid point and we should assign zero for that value
    voxels1 = [x_voxel(counter(j)+1:counter(j+1)),y_voxel(counter(j)+1:counter(j+1)),z_voxel(counter(j)+1:counter(j+1))]; % a group of raw voxels
    BinaryVoxelPortion = voxelInfo(counter(j)+1:counter(j+1))'; % the same group of binarized voxels
    DIST = pdist2(peak_voxels,voxels1); % distance between peak voxels and the group of voxels
    [~,ind] = min(DIST,[],1); % ind stores the index of the assigned nearest peak to each voxel
    ind(BinaryVoxelPortion == 0) = 0; % to assign zero index value to solid voxels
    throat_body(counter(j)+1:counter(j+1)) = ind';
end
% incase some throat gets no voxel, remove throat from pore throat center and neighbour;incase some pore body gets no voxel, remove pore body from locals
i_delete_throats=[];i_delete_locals=[];
for i=1:size(peak_voxels,1)
    if (size(find(throat_body==i),1) == 0 && i<=num_pore_body)
        i_delete_locals=[i_delete_locals;i];
    else if (size(find(throat_body==i),1) == 0 && i>num_pore_body)
            i_delete_throats=[i_delete_throats;i];
        end
    end
end
locals(i_delete_locals,1)=0;
throat_center(i_delete_throats-num_pore_body,1)=0;
neighbour(i_delete_throats-num_pore_body,1)=0;
locals=locals(locals(:,1)~=0,:);
throat_center=throat_center(throat_center(:,1)~=0,:);
neighbour=neighbour(neighbour(:,1)~=0,:);
%----------------------------------------------------------------- END CODE--------------------------------------------------------------------------