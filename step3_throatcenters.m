function [throat_center, neighbour,constriction_length1]  = step3_throatcenters(separator,voxels,ste,locals,num_pore_bodies,neighbour)
% This function finds the individual constriction surfaces, its length and throat center on the constriction surfaces. This is part of Steps 2 and 3 of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
% to find boundary surfaces between pores
indices = unique(nonzeros(separator));
throat_center=[];neighbour=[neighbour zeros(size(neighbour,1),1)];
% a throat center is a voxel that has minimum difference in distance from the adjacent pore body centers
% and the sum of distances to the pore body centers is also minimum,ensuring equidistant and geometrically midway center point at shortest distance
constriction_length1=[];
for i= 1:size(neighbour,1)
    index = neighbour(i,1);j=neighbour(i,2);
    %disp(['Portion ',num2str(i),' out of ',num2str(size(neighbour)), 'processed']);
    locations = find(separator==indices(index));% location of 1 pore voxels
    locations1 = find(separator==indices(j));
    %figure, scatter3(voxels(locations,1),voxels(locations,2),voxels(locations,3),'MarkerEdgeColor',[20/255 98/255 214/255],'MarkerFaceColor',[20/255 98/255 214/255]);hold on;
    %scatter3(voxels(locations1,1),voxels(locations1,2),voxels(locations1,3),'MarkerEdgeColor',[20/255 198/255 214/255],'MarkerFaceColor',[20/255 98/255 214/255])
    kdt_1 = KDTreeSearcher(voxels(locations,:));%build a kdtree of 1st pore volume voxels
    [surfaces_ind1,val1] = knnsearch(kdt_1,voxels(locations1,:),'K',1);%find 2 nearest neighbours(1st one being itself) of each voxel in pore volume 2  in pore voxel1
    surfaces_ind2=unique(locations1(val1<=sqrt(3)*ste));
    surfaces_ind1(val1>sqrt(3)*ste)=[];
    surfaces_ind1=unique(locations(surfaces_ind1));
    surfaces1=voxels(surfaces_ind1,:);
    surfaces2=voxels(surfaces_ind2,:);
    if size(surfaces1,1)==0 && size(surfaces2,1)==0% not neighbours
        neighbour(i,3)=0;
        continue
    end
    if size(surfaces1,1)>0
        %to find the constriction center
        dist_sum1= pdist2(surfaces1,voxels(locals(indices(index)),:))+ pdist2(surfaces1,voxels(locals(indices(j)),:)); % sum of distance between surface voxels and corresponding pore body peaks
        dist_sum1= roundn(dist_sum1,-3);[minDist_sum1,~] = min(dist_sum1,[],1);ind_sum1= find(dist_sum1 == minDist_sum1);
        %centers_plot=[voxels(locals(indices(index)),:); voxels(locals(indices(j)),:); surfaces(ind_sum1,:)]; figure, scatter3(centers_plot(:,1),centers_plot(:,2), centers_plot(:,3))
        dist_diff1 = abs(pdist2(surfaces1(ind_sum1,:),voxels(locals(indices(index)),:)) - pdist2(surfaces1(ind_sum1,:),voxels(locals(indices(j)),:))); % difference of distance between surface voxels with minimum sum and corresponding pore body peaks
        dist_diff1= roundn(dist_diff1,-3);[minDist_diff1,ind_diff1] = min(dist_diff1,[],1);ind_diff3= find(dist_diff1 == minDist_diff1);
        if (size(ind_diff3,1) == 1); center1 = surfaces1(ind_sum1(ind_diff1),:); else; center1 = surfaces1(ind_sum1(ind_diff1(1)),:);end % there shouldnt be more than 2 voxels detected
    end
    % find common surface and center of 2nd pore body surface
    if size(surfaces2,1)>0
        %to find the constriction center
        dist_sum2= pdist2(surfaces2,voxels(locals(indices(index)),:))+ pdist2(surfaces2,voxels(locals(indices(j)),:)); % sum of distance between surface voxels and corresponding pore body peaks
        dist_sum2= roundn(dist_sum2,-3);[minDist_sum2,~] = min(dist_sum2,[],1);ind_sum2= find(dist_sum2 == minDist_sum2);
        %centers_plot=[voxels(locals(indices(index)),:); voxels(locals(indices(j)),:); surfaces(ind_sum1,:)]; figure, scatter3(centers_plot(:,1),centers_plot(:,2), centers_plot(:,3))
        dist_diff2 = abs(pdist2(surfaces2(ind_sum2,:),voxels(locals(indices(index)),:)) - pdist2(surfaces2(ind_sum2,:),voxels(locals(indices(j)),:))); % difference of distance between surface voxels with minimum sum and corresponding pore body peaks
        dist_diff2= roundn(dist_diff2,-3);[minDist_diff2,ind_diff2] = min(dist_diff2,[],1);ind_diff4= find(dist_diff2 == minDist_diff2);
        if (size(ind_diff4,1) == 1); center2 = surfaces2(ind_sum2(ind_diff2),:); else; center2 = surfaces2(ind_sum2(ind_diff2(1)),:);end % there shouldnt be more than 2 voxels detected
    end
    %if both centers coincide anycenter is the center, else the one with minimum difference in dist between 2 body centers is the center, else just pick any center among the 2
    if size(surfaces1,1)>0 && size(surfaces2,1)>0
        if center1(:)==center2(:)
            center=center1;surfaces=surfaces1;surfaces_ind=surfaces_ind1;ind_sum=ind_sum1;ind_diff=ind_diff1;
        else if abs(pdist2(center1,voxels(locals(indices(index)),:)) - pdist2(center1,voxels(locals(indices(j)),:))) < abs(pdist2(center2,voxels(locals(indices(index)),:)) - pdist2(center2,voxels(locals(indices(j)),:)))
                center=center1; surfaces=surfaces1;surfaces_ind=surfaces_ind1;ind_sum=ind_sum1;ind_diff=ind_diff1;
            else if abs(pdist2(center1,voxels(locals(indices(index)),:)) - pdist2(center1,voxels(locals(indices(j)),:))) > abs(pdist2(center2,voxels(locals(indices(index)),:)) - pdist2(center2,voxels(locals(indices(j)),:)))
                    center=center2;surfaces=surfaces2;surfaces_ind=surfaces_ind2;ind_sum=ind_sum2;ind_diff=ind_diff2;
                else
                    center=center1; surfaces=surfaces1;surfaces_ind=surfaces_ind1;ind_sum=ind_sum1;ind_diff=ind_diff1;
                end
            end
        end
    else if size(surfaces1,1)>0
            center=center1;surfaces=surfaces1;surfaces_ind=surfaces_ind1;ind_sum=ind_sum1;ind_diff=ind_diff1;
        else
            center=center2;surfaces=surfaces2;surfaces_ind=surfaces_ind2;ind_sum=ind_sum2;ind_diff=ind_diff2;
        end
    end
    % consider only those throats as throats whose constriction sizes are less than the distance between the adjacent pore bodies
    [~,ind]=max(pdist2(center,surfaces,'euclidean'));
    length1=max(pdist2(surfaces(ind,:),surfaces,'euclidean'));
    if (length1/2 < min(pdist2(center,voxels(locals(indices(index)),:)), pdist2(center,voxels(locals(indices(j)),:)))) && (length1 > 5*ste)
        throat_center = [throat_center; center surfaces_ind(ind_sum(ind_diff)) indices(index) indices(j)];% throat_center, its voxel coordinate, its pore body center coordinates
        %hold on; scatter3(center(1,1),center(1,2),center(1,3),'r','filled');hold on
        %scatter3(surfaces(:,1),surfaces(:,2),surfaces(:,3),'y','filled')
        constriction_length1=[constriction_length1; length1];
        neighbour(i,3)=num_pore_bodies+1;
        num_pore_bodies=num_pore_bodies+1;
    end
    %end
    %centers_plot=[voxels(locals(indices(index)),:); voxels(locals(indices(j)),:); center]; figure, scatter3(centers_plot(:,1),centers_plot(:,2), centers_plot(:,3))
end
%centers_plot=voxels(locals,:); figure, scatter3(centers_plot(:,1),centers_plot(:,2), centers_plot(:,3))
%centers_plot=throat_center(:,1:3); figure, scatter3(centers_plot(:,1),centers_plot(:,2), centers_plot(:,3))
neighbour=neighbour(neighbour(:,3)>0,:);
%----------------------------------------------------------------- END CODE--------------------------------------------------------------------------