function [firstq,secondq,thirdq,fourthq,firstq1,secondq1,thirdq1,fourthq1,throat_axis_center]=step5g_inflexion(store_new,surface_profile_xy1,vis)
% This function finds candidate inflexion points if present along the throat boundary in the 4 quadrants. it is part of Step 5g of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
l=sortrows(knnsearch(store_new(:,1:2),surface_profile_xy1));surface_profile_xy1=[];k=1;l1=[];l=unique(l);
for i=1:size(l,1)
    if i== size(l,1)
        surface_profile_xy1=[surface_profile_xy1;store_new(l(k):l(i),1:2)];l1=[l1;(l(k):l(i)).'];
        break
    end
    if abs(l(i,1)-l(i+1,1))<=round(size(store_new,1)/500) || abs(store_new(l(i),1)-store_new(l(i+1),1))<=min(abs(max(store_new(:,1))-min(store_new(:,1))),abs(max(store_new(:,2))-min(store_new(:,2))))/50% to find nearest point of surface projection to computed boundary based on throat dimension- distance<width of throat/20
        continue
    else
        surface_profile_xy1=[surface_profile_xy1;store_new(l(k):l(i),1:2)];l1=[l1;(l(k):l(i)).'];k=i+1;
    end
end

slope1=[];k=2;xs=store_new(1,1);ys=store_new(1,2);xm=store_new(1+k,1);ym=store_new(1+k,2);% intializations
for p=2:(size(store_new,1)-2)
    if p*k+1<size(store_new,1)
        xe=store_new(p*k+1,1);ye=store_new(p*k+1,2);ys1=(ym-ys)/(xm-xs);
        ys2=(ye-ym)/(xe-xm);y1=(ys1+ys2)/2;yy=atand(y1);h = xe-xm;Nd=xm-xs;Nay=Nd*(ye-ym)+h*(ys-ym);
        Nax=(Nd*h*(Nd+h))/2;y2=Nay/Nax;slope1=[slope1; y2 xs ys xm ym xe ye];%y2 is the second derivative at that point
        xs=xm;ys=ym;xm=xe;ym=ye;
    else
        break;
    end
end
infl=[];

% to find the curvature changing points or inflexion points in the upper 2 parts of the boundary, im finding points that have different sign of slope derivatives
%on either side for a long time
for i = 2: size(slope1,1)-1
    if size(surface_profile_xy1,1)>0 && ismember([slope1(i,2),slope1(i,3)],surface_profile_xy1,'rows')~=1 && ismember([slope1(i,4),slope1(i,5)],surface_profile_xy1,'rows')~=1 && ismember([slope1(i,6),slope1(i,7)],surface_profile_xy1,'rows')~=1
        if i<4
            if sign(slope1(i-1,1))~= sign(slope1(i+1,1)) && sign(slope1(i-1,1))~= sign(slope1(i+2,1)) && sign(slope1(i-1,1))~= sign(slope1(i+3,1)) && abs(slope1(i,1)) < 400 %&& sum(slope1(1:i-1,1)>500)==0 && sum(slope1(end-5+i:end,1)>500)==0
                infl=[infl; slope1(i,:) ];
            end
        else
            if i> size(slope1,1)-4;if sign(slope1(i-1,1))~= sign(slope1(i+1,1)) && sign(slope1(i-2,1))~= sign(slope1(i+1,1)) && sign(slope1(i-3,1))~= sign(slope1(i+1,1))&& abs(slope1(i,1)) < 400  %&& sum(slope1(i-5:i-1,1)>500)==0
                    infl=[infl; slope1(i,:) ];
                end
            else; if sign(slope1(i-1,1))~= sign(slope1(i+1,1)) && sign(slope1(i-2,1))~= sign(slope1(i+1,1)) && sign(slope1(i-3,1))~= sign(slope1(i+3,1)) && abs(slope1(i,1)) < 400  %&& sum(slope1(i-5:i-1,1)>500)==0
                    infl=[infl; slope1(i,:) ];
                end
            end
        end
    else% if we have no surface profile for a certain throat-body combination, then we go ahead and find inflexion points along the entire boundary of throat
        if i<4;if sign(slope1(i-1,1))~= sign(slope1(i+1,1)) && sign(slope1(i-1,1))~= sign(slope1(i+2,1)) && sign(slope1(i-1,1))~= sign(slope1(i+3,1)) && abs(slope1(i,1)) < 400 %&& sum(slope1(1:i-1,1)>500)==0 && sum(slope1(end-5+i:end,1)>500)==0
                infl=[infl; slope1(i,:) ];
            end
        else if i> size(slope1,1)-4;if sign(slope1(i-1,1))~= sign(slope1(i+1,1)) && sign(slope1(i-2,1))~= sign(slope1(i+1,1)) && sign(slope1(i-3,1))~= sign(slope1(i+1,1))&& abs(slope1(i,1)) < 400  %&& sum(slope1(i-5:i-1,1)>500)==0
                    infl=[infl; slope1(i,:) ];
                end
            else; if sign(slope1(i-1,1))~= sign(slope1(i+1,1)) && sign(slope1(i-2,1))~= sign(slope1(i+1,1)) && sign(slope1(i-3,1))~= sign(slope1(i+3,1)) && abs(slope1(i,1)) < 400  %&& sum(slope1(i-5:i-1,1)>500)==0
                    infl=[infl; slope1(i,:) ];
                end
            end
        end
    end
end

% post processing to find the cornermost 4 inflexion points based on first inflexion point going outwards towards body centers taking length of throat as the axis
firstq=[];secondq=[];thirdq=[];fourthq=[];throat_axis_center=[mean(store_new(:,1)),mean(store_new(:,2))];%hold on,plot(throat_axis_center(:,1),throat_axis_center(:,2),'m*')
if size(infl,1)>0;if vis ==1;hold on;plot(infl(:, 2),infl(:, 3), 'm*');hold on;plot(infl(:, 4),infl(:, 5), 'm*');hold on;plot(infl(:, 6),infl(:, 7), 'm*');hold on;end
    infl_coordinates=[infl(:,2:3);infl(:,4:5);infl(:,6:7)];infl_coordinates= unique(infl_coordinates(:,1:2), 'rows', 'stable');%to remove duplicate values from infi1 and find unique 4 values
    %hold on,plot(throat_axis_center(:,1),throat_axis_center(:,2),'m*');
    % to get the inflexion points from the extreme points close to the body centers
    for i=1: size(infl_coordinates,1)% to get the inflexion points from the extreme points close to the body centers
        if infl_coordinates(i,1)>= throat_axis_center(1,1) && infl_coordinates(i,2)>= throat_axis_center(1,2)
            firstq=[firstq; infl_coordinates(i,:)];
        end
        if infl_coordinates(i,1)< throat_axis_center(1,1) && infl_coordinates(i,2)> throat_axis_center(1,2)
            secondq=[secondq; infl_coordinates(i,:)];
        end
        if infl_coordinates(i,1)< throat_axis_center(1,1) && infl_coordinates(i,2)< throat_axis_center(1,2)
            thirdq=[thirdq; infl_coordinates(i,:)];
        end
        if infl_coordinates(i,1)>= throat_axis_center(1,1) && infl_coordinates(i,2)<= throat_axis_center(1,2)
            fourthq=[fourthq; infl_coordinates(i,:)];
        end
    end
end

% incase the throat or combined profile , none have inflexion points in some quardrant , we need to find highest or lowest points of the quadrants
firstq1=[];secondq1=[];thirdq1=[];fourthq1=[];
for i=1: size(store_new,1)% to get the inflexion points from the extreme points close to the body centers
    if store_new(i,1)>= throat_axis_center(1,1) && store_new(i,2)>= throat_axis_center(1,2);firstq1=[firstq1; store_new(i,1:2)];end
    if store_new(i,1)< throat_axis_center(1,1) && store_new(i,2)> throat_axis_center(1,2);secondq1=[secondq1; store_new(i,1:2)];end
    if store_new(i,1)< throat_axis_center(1,1) && store_new(i,2)< throat_axis_center(1,2);thirdq1=[thirdq1; store_new(i,1:2)];end
    if store_new(i,1)>= throat_axis_center(1,1) && store_new(i,2)<= throat_axis_center(1,2);fourthq1=[fourthq1; store_new(i,1:2)];end
end

% to further refine the inflexion points in each quadrant so that in 1st quadrant I get rightmost continous point not on the vertical
% point1=[];point2=[];point3=[];point4=[];
% if size(firstq)>0
%     firstq=sortrows(firstq,1);
%     for i=1:size(firstq,1)
%         if i== size(firstq,1);point1=[firstq(i,:) 0];%scatter(point1(1, 1),point1(1, 2), 'r','filled');hold on;
%             break;end
%         if pdist2(firstq(i+1,:),firstq(i,:))<tol;continue;else;point1=[firstq(i,:) 0];%scatter(point1(1, 1),point1(1, 2), 'r','filled');hold on;
%             break;end;end;end
% if size(secondq)>0
%     secondq=sortrows(secondq,1,'descend');
%     for i=1:size(secondq,1)
%         if i== size(secondq,1);point2=[secondq(i,:) 0];%scatter(point2(1, 1),point2(1, 2), 'r','filled');hold on;
%             break;end
%         if pdist2(secondq(i,:),secondq(i+1,:))<tol;continue;else;point2=[secondq(i,:) 0];%scatter(point2(1, 1),point2(1, 2), 'r','filled');hold on;
%             break;end;end;end
% if size(thirdq)>0
%     thirdq=sortrows(thirdq,1,'descend');
%     for i=1:size(thirdq,1)
%         if i== size(thirdq,1);point3=[thirdq(i,:) 0];%scatter(point3(1, 1),point3(1, 2), 'r','filled');hold on;
%             break;end
%         if thirdq(i,1)-thirdq(i+1,1)<tol;continue;else;point3=[thirdq(i,:) 0];%scatter(point3(1, 1),point3(1, 2), 'r','filled');hold on;
%             break;end;end;end
% if size(fourthq)>0
%     fourthq=sortrows(fourthq,1);
%     for i=1:size(fourthq,1)
%         if i== size(fourthq,1);point4=[fourthq(i,:) 0];%scatter(point4(1, 1),point4(1, 2), 'r','filled');hold on;
%             break;end
%         if fourthq(i+1,1)-fourthq(i,1)<tol;continue;else;point4=[fourthq(i,:) 0];%scatter(point4(1, 1),point4(1, 2), 'r','filled');hold on;
%             break;end;end;end
%--------------------------------------------------------------------------END CODE-----------------------------------------------------------------------------------