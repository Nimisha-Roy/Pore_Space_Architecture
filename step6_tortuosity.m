function [tortuosity,shortest_paths]=step6_tortuosity(throat_body_new,throat_center_new,unconnected_pores,voxels,num_pore_body,entry_direction,locals_new,sTI,neighbour,ste)
% This function computed the tortuosity of all possible paths along the input direction in the structure
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
% to detect all the unconnected pores and throats to form a connected graph
% to find the pore body that has 1 throat (dead end)
if size(unconnected_pores,1)>0
    ind=arrayfun(@(x) find(throat_center_new(:,5)==x,1),unconnected_pores,'un',0);ind=ind(~cellfun('isempty',ind));%find dead end pores in column 5
    ind1=arrayfun(@(x) find(throat_center_new(:,6)==x,1),unconnected_pores,'un',0);ind1=ind1(~cellfun('isempty',ind1));%find dead end pores in column 6
    neighbour_new=[];neighbour_new1=[];
    for a=1:size(throat_center_new,1)% to find distances between neighbouring pore body and throats
        if locals_new(throat_center_new(a,5),1)~=0 && locals_new(throat_center_new(a,6),1)~=0% if size of pore body is 0, locals_new is 0 and so we are assigning 0 distance here and will later remove these from computations
            neighbour_new=[neighbour_new;neighbour(a,3) neighbour(a,1) (pdist2(voxels(throat_center_new(a,4),:),voxels(locals_new(throat_center_new(a,5)),:)))];
            neighbour_new1=[neighbour_new1;neighbour(a,3) neighbour(a,2) (pdist2(voxels(throat_center_new(a,4),:),voxels(locals_new(throat_center_new(a,6)),:)))];
        else
            neighbour_new=[neighbour_new;neighbour(a,3) neighbour(a,1) 0];
            neighbour_new1=[neighbour_new1;a+neighbour(a,3) neighbour(a,2) 0];
        end
    end
    % throat-pore in column6-distance between centers
    if size(cell2mat(ind),1)>0
        neighbour_new(cell2mat(ind),:)=[];% remove the dead end pores in column 5
    end
    if size(cell2mat(ind1),1)>0
        neighbour_new1(cell2mat(ind1),:)=[];% remove the dead end pores in column 6
    end
    neighbour_new=neighbour_new(neighbour_new(:,3)~=0,:);neighbour_new1=neighbour_new1(neighbour_new1(:,3)~=0,:);
    neighbour_pore_throat=[neighbour_new;neighbour_new1];% final neighbour information for input to graph
else
    neighbour_new=[];neighbour_new1=[];
    for a=1:size(throat_center_new,1)% to find distances between neighbouring pore body and throats
        if locals_new(throat_center_new(a,5),1)~=0 && locals_new(throat_center_new(a,6),1)~=0% if size of pore body is 0, locals_new is 0 and so we are assigning 0 distance here and will later remove these from computations
            neighbour_new=[neighbour_new;neighbour(a,3) neighbour(a,1) (pdist2(voxels(throat_center_new(a,4),:),voxels(locals_new(throat_center_new(a,5)),:)))];
            neighbour_new1=[neighbour_new1;neighbour(a,3) neighbour(a,2) (pdist2(voxels(throat_center_new(a,4),:),voxels(locals_new(throat_center_new(a,6)),:)))];
        else
            neighbour_new=[neighbour_new;a+num_pore_body throat_center_new(a,4) 0];
            neighbour_new1=[neighbour_new1;a+num_pore_body throat_center_new(a,5) 0];
        end
    end
    neighbour_pore_throat=[neighbour_new;neighbour_new1];
end
%remove 0 weight edges from neigbours
neighbour_pore_throat(neighbour_pore_throat(:,3)==0,:)=[];
neighbour_pore_throat=unique(neighbour_pore_throat,'rows');

%% graph construction and finding entry and exit surface pores
G = graph(neighbour_pore_throat(:,1),neighbour_pore_throat(:,2),neighbour_pore_throat(:,3)) ;
%figure,plot(G);
if entry_direction==3
    throat_body_new1=throat_body_new;
else if entry_direction==1
        throat_body_new1=permute(throat_body_new,[3,2,1]);
    else
        throat_body_new1=permute(throat_body_new,[1,3,2]);
    end
end
entry_surface_pores=unique(nonzeros(throat_body_new1(1:sTI(1)*sTI(2))));%entry surface pores
exit_surface_pores=unique(nonzeros(throat_body_new1(size(voxels,1)-sTI(1)*sTI(2)+1:size(voxels,1))));%exit surface pores

% remove unconnected pores from entry and exit surface
ind=arrayfun(@(x) find(entry_surface_pores==x),unconnected_pores,'un',0);%find dead end pores in entry surface
ind1=arrayfun(@(x) find(exit_surface_pores==x),unconnected_pores,'un',0);%find dead end pores in exit surface

if size(cell2mat(ind),1)>0
    entry_surface_pores(cell2mat(ind),:)=[];% remove the dead end pores in column 5
end
if size(cell2mat(ind1),1)>0
    exit_surface_pores(cell2mat(ind1),:)=[];% remove the dead end pores in column 6
end
%add equal weight edges to the exit surface pores to obtain the sink node
for i=1:size(exit_surface_pores,1)
    G=addedge(G,exit_surface_pores(i,1),(size(throat_center_new,1)+num_pore_body+1),1);
    %disp(['Portion of ',num2str(i),' out of ',num2str(length(exit_surface_pores)), 'processed']);
end

%% finding shortest path lengths of each entry surface pore to some exit surface pore
tortuosity=[];%figure,p=plot(G);
shortest_paths=zeros(size(entry_surface_pores,1),size(entry_surface_pores,1));count=0;all_paths=[];
%p=figure,
for i=1:size(entry_surface_pores,1)-1
    sp=shortestpath(G,entry_surface_pores(i),(size(throat_center_new,1)+num_pore_body+1));
    shortest_paths(i,1:size(sp,2)-2)=sp(2:end-1);
    if size(sp,1)==0
        continue
    end
    d=distances(G,entry_surface_pores(i),(size(throat_center_new,1)+num_pore_body+1)) ;
    TR = shortestpathtree(G,entry_surface_pores(i),(size(throat_center_new,1)+num_pore_body+1));
    %highlight(p,TR,'EdgeColor','r', 'LineWidth',3);hold on;
    if sp(1)<=num_pore_body%find entry and exit pore center coordinate
        entry_pore_coordinate=voxels(locals_new(sp(1)),:);
    else
        entry_pore_coordinate=voxels(throat_center_new(sp(1)-num_pore_body,4),:);
    end
    if sp(end-1)<=num_pore_body
        exit_pore_coordinate=voxels(locals_new(sp(end-1)),:);
    else
        exit_pore_coordinate=voxels(throat_center_new(sp(end-1)-num_pore_body,4),:);
    end
    path_coord=[];
    for j=1:size(sp,2)-1
        if sp(j)<=num_pore_body%find entry and exit pore center coordinate
            path_coord=[path_coord;voxels(locals_new(sp(j)),:)];
        else
            path_coord=[path_coord;throat_center_new(neighbour(:,3)==sp(j),1:3)];
        end
    end
    all_paths=[all_paths;path_coord/ste];
    count=count+1;
    %plot3(path_coord(:,1),path_coord(:,2),path_coord(:,3),'MarkerSize',5,'MarkerEdgeColor','k','Color',[165/255,1/255,1/255],'Linewidth',2,'MarkerSize',2);hold on
    tortuosity=[tortuosity; (d-1)/pdist2(entry_pore_coordinate,exit_pore_coordinate)];
    %disp(['Portion of ',num2str(i),' out of ',num2str(length(entry_surface_pores)), 'processed']);
end
%ax=gca;
%ax.ZDir='reverse';set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);set(gca,'ZTickLabel',[]);set(gca, 'TickLength',[0 0]);
%writematrix(all_paths,'paths.csv')
%check that no tortuosity value lesser than 1
%sum(tortuosity<1)
%----------------------------------------------------------------- END CODE--------------------------------------------------------------------------

