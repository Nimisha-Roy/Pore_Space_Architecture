function throat_body_new=step5i_correct_division3d(x1,throat_body_new,indices,index,throat_center,body_center_xy1,points1,body1_center,body2_center,center,unit_norm,combined_profile_i)
% This function finds the correct division of pore volume into pore bodies and throats based on 3D inflexion points partitioning. This is part of step 5i of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
ind=find(throat_body_new==indices(index));%to find indicies of current throat voxels
ind_b1=find(throat_body_new==throat_center(x1,5));ind_b2=find(throat_body_new==throat_center(x1,6));% to find indices of adjacent body voxels too
ind_all=[ind;ind_b1;ind_b2];%pt=voxels(locals(throat_center(x1,5)),:);v1=points1(1,:);v4=points1(4,:);v2=points1(2,:);v3=points1(3,:);
[~,ind1]=min(body_center_xy1(:,1));
% find which side of the line connecting positions 1and4 is the
% throatn center by taking sign of cross product and assign the opposite part to respective body centers
if ind1==2% positions 1 and 4 here with bodycenter 1  and position 3 and 2 with body center 2
    % for points 1 and 4
    AB=points1(4,:)-points1(1,:);AC2=body2_center-points1(1,:);
    db2= sign(sind(atan2d(dot(unit_norm.', cross(AB,AC2)),dot(AB,AC2))));% find which side of the line connecting positions 1and4 is throat center,bodycenter 1 and 2 by taking sign of cross product
    ACt=center-points1(1,:);
    dbt= sign(sind(atan2d(dot(unit_norm.', cross(AB,ACt)),dot(AB,ACt))));
    AC1=body1_center-points1(1,:);
    db1= sign(sind(atan2d(dot(unit_norm.', cross(AB,AC1)),dot(AB,AC1))));
    AB1=(repmat(points1(4,:),size(ind_all,1),1)-repmat(points1(1,:),size(ind_all,1),1));
    AP= (combined_profile_i-repmat(points1(1,:),size(ind_all,1),1));
    d1= sign(sind(atan2d(dot(repmat(unit_norm.',size(ind_all,1),1),cross(AB1,AP),2),dot(AB1,AP,2))));% sign of all voxels with respect to points 1 and points 4
    % for points 2 and 3
    ABi=points1(3,:)-points1(2,:);AC2i=body2_center-points1(1,:);
    db2i= sign(sind(atan2d(dot(unit_norm.', cross(ABi,AC2i)),dot(ABi,AC2i))));% find which side of the line connecting positions 2and3 is throat center,bodycenter 1 and 2 by taking sign of cross product
    ACti=center-points1(2,:);
    dbti= sign(sind(atan2d(dot(unit_norm.', cross(ABi,ACti)),dot(ABi,ACti))));
    AC1i=body1_center-points1(2,:);
    db1i= sign(sind(atan2d(dot(unit_norm.', cross(ABi,AC1i)),dot(ABi,AC1i))));
    AB1i=(repmat(points1(3,:),size(ind_all,1),1)-repmat(points1(2,:),size(ind_all,1),1));
    APi= (combined_profile_i-repmat(points1(2,:),size(ind_all,1),1));
    d1i= sign(sind(atan2d(dot(repmat(unit_norm.',size(ind_all,1),1),cross(AB1i,APi),2),dot(AB1i,APi,2))));% sign of all voxels with respect to points 2 and points 3
    
    if sign(db1)~=sign(db2i) && (sign(dbt)~=sign(db1)|| sign(dbti)~=sign(db2i)) % includes 3 case scenarios .|.|. .| |.. ..||.
        throat_body_new(ind_all(d1(:,1)==db1))=throat_center(x1,6);%assign all voxels on right of 1&4 line as body center2
        throat_body_new(ind_all(d1i(:,1)==db2i))=throat_center(x1,5);%assign all voxels on left of 2&3 line as body center1
        throat_body_new(ind_all(d1(:,1)~=db1 & d1i(:,1)~=db2i))=indices(index);%assign all voxels on right of 2&3 and left of 1&4 as throat
        
    else
        if sign(db1)~=sign(db2i) && (sign(dbt)==sign(db1)&& sign(dbti)==sign(db2i))% includes 1 case scenario |...|
            throat_body_new(ind_all(d1(:,1)~=db1))=throat_center(x1,6);
            throat_body_new(ind_all(d1i(:,1)~=db2i))=throat_center(x1,5);
            throat_body_new(ind_all(d1(:,1)==db1 & d1i(:,1)==db2i))=indices(index);
        else
            if sign(db1)==sign(db2i)
                
                if sign(dbti)==sign(db2i) && sign(dbt)~=sign(db1) % case |..|.
                    throat_body_new(ind_all(d1(:,1)==db1))=throat_center(x1,6);
                    throat_body_new(ind_all(d1i(:,1)~=db2i))=throat_center(x1,5);
                    throat_body_new(ind_all(d1(:,1)~=db1 & d1i(:,1)==db2i))=indices(index);
                else
                    if sign(dbti)~=sign(db2i) && sign(dbt)==sign(db1)% case .|..|
                        throat_body_new(ind_all(d1(:,1)~=db1))=throat_center(x1,6);
                        throat_body_new(ind_all(d1i(:,1)==db2i))=throat_center(x1,5);
                        throat_body_new(ind_all(d1(:,1)==db1 & d1i(:,1)~=db2i))=indices(index);
                    else
                        if sign(dbti)==sign(db2i) && sign(dbt)==sign(db1) && sign(dbt)~= sign(db2)% case |.|..
                            throat_body_new(ind_all(d1(:,1)==db1))=throat_center(x1,6);
                            throat_body_new(ind_all(d1i(:,1)~=db2i))=throat_center(x1,5);
                            throat_body_new(ind_all(d1(:,1)~=db1 & d1i(:,1)==db2i))=indices(index);
                        else
                            if sign(dbti)==sign(db2i) && sign(dbt)==sign(db1) && sign (dbt)==sign(db2)% case ..|.|
                                throat_body_new(ind_all(d1(:,1)~=db1))=throat_center(x1,6);
                                throat_body_new(ind_all(d1i(:,1)==db2i))=throat_center(x1,5);
                                throat_body_new(ind_all(d1(:,1)==db1 & d1i(:,1)~=db2i))=indices(index);
                            end
                        end
                    end
                end
            end
        end
    end
else% positions 3 and 2 here with body center 1 and 1 and 4 with 2
    % for points 1 and 4
    AB=points1(4,:)-points1(1,:);
    AC2=body2_center-points1(1,:);
    db2= sign(sind(atan2d(dot(unit_norm.', cross(AB,AC2)),dot(AB,AC2))));% find which side of the line connecting positions 1and4 is throat center,bodycenter 1 and 2 by taking sign of cross product
    ACt=center-points1(1,:);
    dbt= sign(sind(atan2d(dot(unit_norm.', cross(AB,ACt)),dot(AB,ACt))));
    AC1=body1_center-points1(1,:);
    db1= sign(sind(atan2d(dot(unit_norm.', cross(AB,AC1)),dot(AB,AC1))));
    AB1=(repmat(points1(4,:),size(ind_all,1),1)-repmat(points1(1,:),size(ind_all,1),1));
    AP= (combined_profile_i-repmat(points1(1,:),size(ind_all,1),1));
    d1= sign(sind(atan2d(dot(repmat(unit_norm.',size(ind_all,1),1),cross(AB1,AP),2),dot(AB1,AP,2))));% sign of all voxels with respect to points 1 and points 4
    % for points 2 and 3
    ABi=points1(3,:)-points1(2,:);
    AC2i=body2_center-points1(2,:);
    db2i= sign(sind(atan2d(dot(unit_norm.', cross(ABi,AC2i)),dot(ABi,AC2i))));% find which side of the line connecting positions 2and3 is throat center,bodycenter 1 and 2 by taking sign of cross product
    ACti=center-points1(2,:);
    dbti= sign(sind(atan2d(dot(unit_norm.', cross(ABi,ACti)),dot(ABi,ACti))));
    AC1i=body1_center-points1(2,:);
    db1i= sign(sind(atan2d(dot(unit_norm.', cross(ABi,AC1i)),dot(ABi,AC1i))));
    AB1i=(repmat(points1(3,:),size(ind_all,1),1)-repmat(points1(2,:),size(ind_all,1),1));
    APi= (combined_profile_i-repmat(points1(2,:),size(ind_all,1),1));
    d1i= sign(sind(atan2d(dot(repmat(unit_norm.',size(ind_all,1),1),cross(AB1i,APi),2),dot(AB1i,APi,2))));% sign of all voxels with respect to points 2 and points 3
    if sign(db2)~=sign(db1i) && (sign(dbt)~=sign(db2)|| sign(dbti)~=sign(db1i)) % includes 3 case scenarios .|.|. .| |.. ..||.
        throat_body_new(ind_all(d1(:,1)==db2))=throat_center(x1,6);%assign all voxels on right of 1&4 line as body center2
        throat_body_new(ind_all(d1i(:,1)==db1i))=throat_center(x1,5);%assign all voxels on left of 2&3 line as body center1
        throat_body_new(ind_all(d1(:,1)~=db2 & d1i(:,1)~=db1i))=indices(index);%assign all voxels on right of 2&3 and left of 1&4 as throat
        
    else
        if sign(db2)~=sign(db1i) && (sign(dbt)==sign(db2)&& sign(dbti)==sign(db1i))% includes 1 case scenario |...|
            throat_body_new(ind_all(d1(:,1)~=db2))=throat_center(x1,6);
            throat_body_new(ind_all(d1i(:,1)~=db1i))=throat_center(x1,5);
            throat_body_new(ind_all(d1(:,1)==db2 & d1i(:,1)==db1i))=indices(index);
        else
            if sign(db2)==sign(db1i)
                if sign(dbti)==sign(db1i) && sign(dbt)~=sign(db2) % case |..|.
                    throat_body_new(ind_all(d1(:,1)==db2))=throat_center(x1,6);
                    throat_body_new(ind_all(d1i(:,1)~=db1i))=throat_center(x1,5);
                    throat_body_new(ind_all(d1(:,1)~=db2 & d1i(:,1)==db1i))=indices(index);
                else
                    if sign(dbti)~=sign(db1i) && sign(dbt)==sign(db2)% case .|..|
                        throat_body_new(ind_all(d1(:,1)~=db2))=throat_center(x1,6);
                        throat_body_new(ind_all(d1i(:,1)==db1i))=throat_center(x1,5);
                        throat_body_new(ind_all(d1(:,1)==db2 & d1i(:,1)~=db1i))=indices(index);
                    else
                        if sign(dbti)==sign(db1i) && sign(dbt)==sign(db2) && sign(dbt)~= sign(db1)% case |.|..
                            throat_body_new(ind_all(d1(:,1)==db2))=throat_center(x1,6);
                            throat_body_new(ind_all(d1i(:,1)~=db1i))=throat_center(x1,5);
                            throat_body_new(ind_all(d1(:,1)~=db2 & d1i(:,1)==db1i))=indices(index);
                        else
                            if sign(dbti)==sign(db1i) && sign(dbt)==sign(db2) && sign (dbt)==sign(db1)% case ..|.|
                                throat_body_new(ind_all(d1(:,1)~=db2))=throat_center(x1,6);
                                throat_body_new(ind_all(d1i(:,1)==db1i))=throat_center(x1,5);
                                throat_body_new(ind_all(d1(:,1)==db2 & d1i(:,1)~=db1i))=indices(index);
                            end
                        end
                    end
                end
            end
        end
    end
end
%---------------------------------------------------------------------------END CODE--------------------------------------------------------------------------