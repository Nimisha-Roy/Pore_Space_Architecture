function [firstq2,secondq2,thirdq2,fourthq2,firstqc,secondqc,thirdqc,fourthqc]=step5g_inflexion_combined(center,store_new1,throat_center_xy1,vis)
% This function finds the candidate inflexion points along the combined throat and body boundary. It is part of Step 5h of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
%figure, plot(store_new1(:,1),store_new1(:,2),'g');hold on
slope1=[];k=2;
xs=store_new1(1,1);ys=store_new1(1,2);xm=store_new1(1+k,1);ym=store_new1(1+k,2);% intializations
for p=2:(size(store_new1,1)-2)
    if p*k+1<size(store_new1,1)
        xe=store_new1(p*k+1,1);ye=store_new1(p*k+1,2);ys1=(ym-ys)/(xm-xs);
        ys2=(ye-ym)/(xe-xm);y1=(ys1+ys2)/2;yy=atand(y1);h = xe-xm;Nd=xm-xs;Nay=Nd*(ye-ym)+h*(ys-ym);
        Nax=(Nd*h*(Nd+h))/2;y2=Nay/Nax;
        slope1=[slope1; y2 xs ys xm ym xe ye];%y2 is the second derivative at that point
        xs=xm;ys=ym;xm=xe;ym=ye;
    else
        break;
    end
end
infl=[];

slope1(:,1)=(round(slope1(:,1)*10^(2)))/10^(2);
i=2;
while(i<=size(slope1,1)-1)
    if i<4
        if sign(slope1(i-1,1))~= sign(slope1(i+1,1)) && sign(slope1(i-1,1))~= sign(slope1(i+2,1)) && sign(slope1(i-1,1))~= sign(slope1(i+3,1)) && abs(slope1(i,1)) < 400 %&& sum(slope1(1:i-1,1)>500)==0 && sum(slope1(end-5+i:end,1)>500)==0
            infl=[infl; slope1(i,:) ]; i=i+4;
            continue;
        end
    else if i> size(slope1,1)-4
            if sign(slope1(i-1,1))~= sign(slope1(i+1,1)) && sign(slope1(i-2,1))~= sign(slope1(i+1,1)) && sign(slope1(i-3,1))~= sign(slope1(i+1,1))&& abs(slope1(i,1)) < 400  %&& sum(slope1(i-5:i-1,1)>500)==0
                infl=[infl; slope1(i,:) ];
                break
            end
        else
            if sign(slope1(i-1,1))~= sign(slope1(i+1,1)) && sign(slope1(i-2,1))~= sign(slope1(i+1,1)) && sign(slope1(i-3,1))~= sign(slope1(i+3,1)) && abs(slope1(i,1)) < 400  %&& sum(slope1(i-5:i-1,1)>500)==0
                infl=[infl; slope1(i,:) ];i=i+4;
                continue
            end
        end
    end
    i=i+1;
end

% post processing to find the cornermost 4 inflexion points based on first inflexion point going outwards towards body centers taking length of throat as the axis
firstq2=[];secondq2=[];thirdq2=[];fourthq2=[];throat_axis_center=center;%hold on,plot(throat_axis_center(:,1),throat_axis_center(:,2),'b*')
firstqc=[];secondqc=[];thirdqc=[];fourthqc=[];
if size(infl,1)>0
    if vis==1
        hold on;plot(infl(:, 2),infl(:, 3), 'c*');hold on;plot(infl(:, 4),infl(:, 5), 'c*');hold on;plot(infl(:, 6),infl(:, 7), 'c*');
    end
    infl_coordinates=cell2mat(arrayfun(@(i) horzcat(reshape(infl(i,2:end),[2,3]).',repmat(infl(i,1),[3,1])),linspace(1,size(infl,1),size(infl,1)),'UniformOutput',false).');
    infl_coordinates= unique(infl_coordinates(:,1:3), 'rows', 'stable');%to remove duplicate values from infi1 and find unique 4 values
    if inpolygon(throat_axis_center(1,1),throat_axis_center(1,2),store_new1(:,1),store_new1(:,2))~=1
        throat_axis_center=throat_center_xy1;%hold on,plot(throat_axis_center(:,1),throat_axis_center(:,2),'b*');
    end
    %thresh=sizy/sizx*0.5;
    %remove flat inflexion points by using dominant points
    x=step5g_inflexion_dominant(infl_coordinates(:,1:2),vis);
    infl_coordinates=infl_coordinates(ismember(infl_coordinates(:,1:2),x,'rows'),:);
    % to get the inflexion points in 4 quadrants
    for i=1: size(infl_coordinates,1)
        if infl_coordinates(i,1)>= throat_axis_center(1,1) && infl_coordinates(i,2)>= throat_axis_center(1,2)
            firstq2=[firstq2; infl_coordinates(i,:)];
        end
        if infl_coordinates(i,1)< throat_axis_center(1,1) && infl_coordinates(i,2)> throat_axis_center(1,2)
            secondq2=[secondq2; infl_coordinates(i,:)];
        end
        if infl_coordinates(i,1)< throat_axis_center(1,1) && infl_coordinates(i,2)< throat_axis_center(1,2)
            thirdq2=[thirdq2; infl_coordinates(i,:)];
        end
        if infl_coordinates(i,1)>= throat_axis_center(1,1) && infl_coordinates(i,2)<= throat_axis_center(1,2)
            fourthq2=[fourthq2; infl_coordinates(i,:)];
        end
    end
    
    %extract combined profile in 4 quadrants to get the protrusion ratio
    for i =1:size(store_new1,1)
        if store_new1(i,1)>= throat_axis_center(1,1) && store_new1(i,2)>= throat_axis_center(1,2)
            firstqc=[firstqc; store_new1(i,1:2)];
        end
        if store_new1(i,1)< throat_axis_center(1,1) && store_new1(i,2)> throat_axis_center(1,2)
            secondqc=[secondqc; store_new1(i,1:2)];
        end
        if store_new1(i,1)< throat_axis_center(1,1) && store_new1(i,2)< throat_axis_center(1,2)
            thirdqc=[thirdqc; store_new1(i,1:2)];
        end
        if store_new1(i,1)>= throat_axis_center(1,1) && store_new1(i,2)<= throat_axis_center(1,2)
            fourthqc=[fourthqc; store_new1(i,1:2)];
        end
    end
    %sort the quadrants in ascending/descing order of x
    if size(firstqc,1)>0
        firstqc=sortrows(firstqc,1);
    end
    if size(fourthqc,1)>0
        fourthqc=sortrows(fourthqc,1);
    end
    if size(secondqc,1)>0
        secondqc=sortrows(secondqc,1,'descend');
    end
    if size(thirdqc,1)>0
        thirdqc=sortrows(thirdqc,1,'descend');
    end
end
%-----------------------------------------------------------------------END CODE-------------------------------------------------------------------------------------