function [throat_body_surface]=step4_throat_body_surface(throat_body,voxels,neighbour)
% This function computes the interface surface voxels between throats and bodies. It is used to distinguish throat-body and throat-particle voxels from the throat boundary. This is a preprocessing step for Step 5 of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
throat_body_surface=[];% a throat body surfaces comprises all voxels that separate throats from their adjacent bodies
for i= 1: size(neighbour,1)
    % for surface between a throat and 1 of its adjacent pore body
    %if throat_center(i,7)==0
    index = neighbour(i,1);j=neighbour(i,3);
    locations = find(throat_body==index);% location of 1 pore voxels
    locations1 = find(throat_body==j);
    kdt_1 = KDTreeSearcher(voxels(locations,:));%build a kdtree of ith pore body voxels
    [surfaces_ind1,val1] = knnsearch(kdt_1,voxels(locations1,:),'K',1);%find 1 nearest neighbour of each voxel in pore body i  in throat voxels i
    surfaces_ind2=unique(locations1(val1<=sqrt(3)));
    surfaces_ind1(val1>sqrt(3))=[];
    surfaces_ind1=unique(locations(surfaces_ind1));
    surfaces=voxels(surfaces_ind1,:);
    surfaces=[surfaces;voxels(surfaces_ind2,:)];
    %figure,scatter3(voxels(locations,1),voxels(locations,2),voxels(locations,3),'b');hold on;scatter3(voxels(locations1,1),voxels(locations1,2),voxels(locations1,3),'y');hold on
    if size(surfaces,1)>0
        throat_index=repmat(neighbour(i,3), size(surfaces,1),1);% repeat value size(surfaces,1) times
        throat_body_surface=[throat_body_surface; surfaces throat_index];% surface, its associated throat index
        surfaces=[];
    end
    % for the surface between the same throat but the other adjacent pore body
    index = neighbour(i,2);
    j=neighbour(i,3);
    locations = find(throat_body==index);% location of 1 pore voxels
    locations1 = find(throat_body==j);
    kdt_1 = KDTreeSearcher(voxels(locations,:));%build a kdtree of 1st pore volume voxels
    [surfaces_ind1,val1] = knnsearch(kdt_1,voxels(locations1,:),'K',1);%find 2 nearest neighbours(1st one being itself) of each voxel in pore volume 2  in pore voxel1
    surfaces_ind2=unique(locations1(val1<=sqrt(3)));
    surfaces_ind1(val1>sqrt(3))=[];
    surfaces_ind1=unique(locations(surfaces_ind1));
    surfaces=voxels(surfaces_ind1,:);
    surfaces=[surfaces;voxels(surfaces_ind2,:)];
    %scatter3(voxels(locations,1),voxels(locations,2),voxels(locations,3),'g');
  if size(surfaces,1)>0
        throat_index=repmat(neighbour(i,3), size(surfaces,1),1);% repeat value size(surfaces,1) times
        throat_body_surface=[throat_body_surface; surfaces throat_index];% surface, its associated throat index
        surfaces=[];
    end
        %disp(['Portion ',num2str(i),' out of ',num2str(size(neighbour)), 'processed']);
        surfaces=[];
end
%----------------------------------------------------------------- END CODE--------------------------------------------------------------------------