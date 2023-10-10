function [throat_body_new,lsi,throat_wall_normalized_roughness,throat_center_new,locals_new]= step5_main(throat_body,voxels,num_pore_body,throat_center,ste,locals,throat_body_surface,sig,sig1,vis)
%This function uses curvature of pore throats and bodies to modify their lengths. This is Step 5 of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------ BEGIN CODE --------------------------------------------------------------------------------------
indices = unique(nonzeros(throat_body));digits(4)%to set precision to 4 decimal places
lsi=[];throat_wall_normalized_roughness=[];
throat_body_new=throat_body;
throat_center_new=zeros(size(throat_center,1),8);throat_center_new(:,1:6)=throat_center;
% a throat mid profile is the projected throat points on the mid best fit plane to it
x1=1;
for index = num_pore_body+1:size(indices,1)% from num_pore_bodies+1 are the stored throats
    % to get projection of throat, body voxels and throat_body_surface on best fit plane in 2D
    locations = find(throat_body_new==indices(index));% location of 1 pore voxel
    if size(locations,1)~=0 % if no throat voxels remain in throat_body_new
        locations_body1=find(throat_body_new==throat_center(x1,5));locations_body2=find(throat_body_new==throat_center(x1,6));
        % printing each throat and its corresponding pore centers
        set_throat=voxels(locations,:);set2=throat_body_surface((throat_body_surface(:,4)==indices(index)),1:3);set2=set2(ismember(set2,set_throat,'rows'),:);
        set_body1=voxels(locations_body1,:);set_body2=voxels(locations_body2,:);
        if vis == 1
            figure, scatter3(set_throat(:,1),set_throat(:,2), set_throat(:,3),'MarkerEdgeColor',[165/255 1/255 1/255],'MarkerFaceColor',[165/255 1/255 1/255]); hold on;%scatter3(set2(:,1),set2(:,2), set2(:,3),'y'); hold on
            scatter3(voxels(locations_body1,1),voxels(locations_body1,2),voxels(locations_body1,3),'MarkerEdgeColor',[20/255 98/255 214/255],'MarkerFaceColor',[20/255 98/255 214/255]);hold on;
            scatter3(voxels(locations_body2,1),voxels(locations_body2,2),voxels(locations_body2,3),'MarkerEdgeColor',[20/255 98/255 214/255],'MarkerFaceColor',[20/255 98/255 214/255]);hold on
            scatter3(voxels(locals(throat_center(x1,5)),1),voxels(locals(throat_center(x1,5)),2), voxels(locals(throat_center(x1,5)),3),'r', 'filled');hold on
            scatter3(voxels(locals(throat_center(x1,6)),1),voxels(locals(throat_center(x1,6)),2), voxels(locals(throat_center(x1,6)),3), 'r', 'filled');hold on
            title ('throat and adjacent pore body voxels in 3D')
        end
        %Finding the best fit plane approximating throat and body voxels and passing through the two body centers in their direction
        % This is the Solution of constrained least squares problem solved using lagrangian function with 2 lagrangian multipliers
        [mid_profile,surface_profile,unit_norm,body1_profile,body2_profile]=step5b_bestfit(set_body1,set_body2,set2,ste,x1,set_throat,voxels,locations,locals,throat_center,vis);
        
        %Projection of throat, adjacent bodies and throat_body_surface from best fit plane on xy, yz or xz plane to do computations in 2D with reference axes
        [mid_profile_xy,body_center_xy,throat_center_xy,surface_profile_xy,body1_profile_xy,body2_profile_xy]=step5c_bestfit2D(locals,throat_center,unit_norm,x1,voxels,mid_profile,surface_profile,body1_profile,body2_profile,vis);
        
        % to rotate all components such that body centers are horizontal
        [mid_profile_xy1,body_center_xy1,throat_center_xy1,surface_profile_xy1,body1_profile_xy1,body2_profile_xy1,angle_rotation]=step5c_rotate2D(sig1,mid_profile_xy,body_center_xy,throat_center_xy,surface_profile_xy,body1_profile_xy,body2_profile_xy,vis);
        %% Extracting boundary of projected throat/combined plot and finding inflexion points along the filtered boundaries
        % considering computation on the current throat only if it has atleast 10 voxels
        if size(mid_profile_xy1,1)>10
            % extract boundary from image generated from plot
            %[boundary,bw,graindata]=plot2image_bound(mid_profile_xy1);
            
            % to extract boundary points at constant degree intervals in polar coordinates in a single layer-taking care of reentrant boundary points
            [store_new]=step5d_boundary(mid_profile_xy1,sig,throat_center_xy1,ste,vis);
            
            % filtering meso and micro features from boundary and getting throat wall roughness
            % From evaluation on different angular and rounded particles, a cutoff of 1.2% of maximum intensity in frequency domain was found consistent with removing
            % roughness features from boundary, doing same analysis to see upward rise,we found a cutoff of 15% for mesolevel cutoff, that cutoff was used to remove meso and micro level features from intensities of the walls of pore throat
            thint=(round(pi/360*10^(sig)))/10^(sig);% 0.5 degrees
            [store_new,throat_wall_normalized_roughness,throat_center_new]=step5ef_filter(x1,store_new,thint,throat_wall_normalized_roughness,throat_center_new,throat_center_xy1,vis);
            
            % find candidate inflexion points if present along the throat boundary in the 4 quadrants
            [firstqti,secondqti,thirdqti,fourthqti,firstqtall,secondqtall,thirdqtall,fourthqtall,throat_axis_center]=step5g_inflexion(store_new,surface_profile_xy1,vis);
            
            %find filtered boundary of the combined profile
            combined_profile= [mid_profile_xy1;body1_profile_xy1;body2_profile_xy1];
            [store_new_combined]= step5ef_filtered_combined(combined_profile,ste,sig,thint,vis);
            if vis==1
                figure,plot(store_new(:,1),store_new(:,2),'g','Linewidth',3);hold on
                scatter(throat_axis_center(1,1),throat_axis_center(1,2),'r','filled');hold on
                if size(firstqti,1)>0;scatter(firstqti(:,1),firstqti(:,2),'r*');end;hold on
                if size(secondqti,1)>0;scatter(secondqti(:,1),secondqti(:,2),'r*');end;hold on
                if size(thirdqti,1)>0;scatter(thirdqti(:,1),thirdqti(:,2),'r*');end;hold on
                if size(fourthqti,1)>0;scatter(fourthqti(:,1),fourthqti(:,2),'r*');end;hold on
                plot(store_new_combined(:,1),store_new_combined(:,2),'b','Linewidth',3);hold on
                title('throat and combined boundary with inflexion points')
            end
            
            % find 4 inflexion points if present along the combined boundary in the 4 quadrants
            sizy=max(store_new(:,2))-min(store_new(:,2));sizx=max(store_new(:,1))-min(store_new(:,1));
            [firstqci,secondqci,thirdqci,fourthqci,firstqcall,secondqcall,thirdqcall,fourthqcall]=step5g_inflexion_combined(throat_axis_center,store_new_combined,throat_center_xy1,vis);
            
            % find the right inflexion points in the 4 quadrants
            %qci and qti are inflexion points of combined profile and throat profile in all quadrants respectfully
            %qcall and qtall are combined boundary points and throat boundary points in 4 quadrants respectively
            [point1,point2,point3,point4]=step5gh_correct_inflexion(sizy,sizx,firstqci,secondqci,thirdqci,fourthqci,firstqtall,secondqtall,thirdqtall,fourthqtall,firstqti,secondqti,thirdqti,fourthqti,firstqcall,secondqcall,thirdqcall,fourthqcall);
            points=[point1;point2;point3;point4];
            if vis ==1
                scatter(points(:,1),points(:,2),'r','filled');axis equal
            end
            %% finding approximate length of throat,normalized roughness & longitudinal sharpness index
            L= max(max(pdist2(points(:,1:2),points(:,1:2),'euclidean')));
            throat_wall_normalized_roughness(end,1)=throat_wall_normalized_roughness(end,1)*100/L;
            throat_center_new(x1,7)=throat_center_new(x1,7)*100/L;
            
            % LSI- if point# has a 0 in its 3rd column, it comes from the throat boundary, else it comes from the combined boundary
            [ang1,ang2]=step5h_lsi(store_new,store_new_combined,points,vis);
            lsi=[lsi;0.5*(ang1+ang2)];
            throat_center_new(x1,8)=0.5*(ang1+ang2);
            
            %% correct division of pore throat and body in 3D
            % finding corresponding 4 inflexion points among pore/throat voxels - 1 possible rotation, 2 back projections, 1 to throat plane and 1 back to throat voxels rotation
            [points1,combined_profile_i]=step5i_inflexion3d(angle_rotation,throat_center_xy1,points,mid_profile_xy1,mid_profile,body1_profile,body2_profile,mid_profile_xy,body1_profile_xy,body2_profile_xy,sig);
            
            %division based on inflexion points in 3d
            center=[mean(combined_profile_i(:,1)),mean(combined_profile_i(:,2)),mean(combined_profile_i(:,3))];
            body1_center=mean(body1_profile(:,1:3));body2_center=mean(body2_profile(:,1:3));
            throat_body_new=step5i_correct_division3d(x1,throat_body_new,indices,index,throat_center,body_center_xy1,points1,body1_center,body2_center,center,unit_norm,combined_profile_i);
            
            %check
            if vis==1
                locations = find(throat_body_new==indices(index));
                locations_body1=find(throat_body_new==throat_center(x1,5));locations_body2=find(throat_body_new==throat_center(x1,6));
                figure,scatter3(voxels(locations,1),voxels(locations,2), voxels(locations,3),'MarkerEdgeColor',[165/255 1/255 1/255],'MarkerFaceColor',[165/255 1/255 1/255]); hold on
                scatter3(voxels(locations_body1,1),voxels(locations_body1,2),voxels(locations_body1,3),'MarkerEdgeColor',[20/255 98/255 214/255],'MarkerFaceColor',[20/255 98/255 214/255]);hold on;
                scatter3(voxels(locations_body2,1),voxels(locations_body2,2),voxels(locations_body2,3),'MarkerEdgeColor',[20/255 98/255 214/255],'MarkerFaceColor',[20/255 98/255 214/255]);hold on
                scatter3(voxels(locals(throat_center(x1,6)),1),voxels(locals(throat_center(x1,6)),2), voxels(locals(throat_center(x1,6)),3), 'r','filled');hold on
                scatter3(voxels(locals(throat_center(x1,5)),1),voxels(locals(throat_center(x1,5)),2), voxels(locals(throat_center(x1,5)),3), 'r','filled');hold on
                title('throat after geometric modification');
            end
            %% update throat and body centers
            ind=find(throat_body_new==indices(index));
            if size(ind,1) ~= 0
                throat_center_new(x1,1:3)=mean(voxels(ind,1:3));
                %compute Euclidean distances of throat center to all center voxels
                B=voxels(ind,1:3);
                distances = sqrt(sum(bsxfun(@minus, B, mean(voxels(ind,1:3))).^2,2));
                throat_center_new(x1,4) = ind(find(distances==min(distances),1,'first'));
                x1=x1+1;
            else
                x1=x1+1;
            end
        else
            x1=x1+1;
        end
    else
        x1=x1+1;
    end
end

%% update final pore body centers
locals_new=[];indices = unique(nonzeros(throat_body));
for index = 1:num_pore_body
    locations_body = find(throat_body_new==indices(index));% location of 1 pore throat voxels
    if size(locations_body,1)>0
        B=voxels(locations_body,1:3);%for body center 1
        distances = sqrt(sum(bsxfun(@minus, B, mean(voxels(locations_body,1:3))).^2,2));
        locals_new=[locals_new;locations_body(find(distances==min(distances),1,'first'))];
    else
        locals_new=[locals_new;0];
    end
end
%---------------------------------------------------------------------- END CODE ----------------------------------------------------------------------------------