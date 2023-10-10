function [store_new]=step5d_boundary(mid_profile_xy1,sig,throat_center_xy,ste,vis)
% This function finds polar coordinates of a closed boundary in a single layer-taking care of reentrant boundary points. It is Step 5d of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
mid_profile_xy1=(round(mid_profile_xy1*10^(sig)))/10^(sig);mid_profile_xy1=unique(mid_profile_xy1,'rows');
bound = boundary(mid_profile_xy1(:,1),mid_profile_xy1(:,2),1);% approximate boundary point coordinates to emsure boundary coordinates plus a bunch of infiltrated coordinates
bound= unique(bound);new_bound=[mid_profile_xy1(bound,1),mid_profile_xy1(bound,2)];

if vis ==1
    %figure,scatter(mid_profile_xy1(:,1),mid_profile_xy1(:,2),'r');hold on
    figure,scatter(new_bound(:,1),new_bound(:,2),'y','filled');hold on
end

% to filter out infiltrated coordinates, ensuring only a layer of boundary coordinates
firstsecondq=[];thirdfourthq=[];
for i=1: size(new_bound,1)% to get the inflexion points from the extreme points close to the body centers
    if (new_bound(i,1)>= throat_center_xy(1,1) && new_bound(i,2)>= throat_center_xy(1,2))||(new_bound(i,1)< throat_center_xy(1,1) && new_bound(i,2)> throat_center_xy(1,2))
        firstsecondq=[firstsecondq; new_bound(i,:)];
    end
    if (new_bound(i,1)< throat_center_xy(1,1) && new_bound(i,2)< throat_center_xy(1,2))||(new_bound(i,1)>= throat_center_xy(1,1) && new_bound(i,2)<= throat_center_xy(1,2))
        thirdfourthq=[thirdfourthq; new_bound(i,:)];
    end
end
firstsecondxmax= max(firstsecondq(:,1));firstsecondxmin= min(firstsecondq(:,1));thirdfourthxmax= max(thirdfourthq(:,1));thirdfourthxmin= min(thirdfourthq(:,1));

% points with same x but different close ys (close ys so that we can deal with reentrant angle effects) are merged into 1 point with highest y in top part and lowest part in bottom part
ind1 = (abs((round(firstsecondq(:,1)*10^(2))/10^(2))-(round(firstsecondxmax(:,1)*10^(2))/10^(2)))>ste);new_bound=firstsecondq(~ind1,:);firstsecondq=firstsecondq(ind1,:);
ind2 = (abs((round(firstsecondq(:,1)*10^(2))/10^(2))-(round(firstsecondxmin(:,1)*10^(2))/10^(2)))>ste);new_bound=[new_bound;firstsecondq(~ind2,:)];firstsecondq=firstsecondq(ind2,:);
ind1 = (abs((round(thirdfourthq(:,1)*10^(2))/10^(2))-(round(thirdfourthxmax(:,1)*10^(2))/10^(2)))>ste);new_bound=[new_bound;thirdfourthq(~ind1,:)];thirdfourthq=thirdfourthq(ind1,:);
ind2 = (abs((round(thirdfourthq(:,1)*10^(2))/10^(2))-(round(thirdfourthxmin(:,1)*10^(2))/10^(2)))>ste);new_bound=[new_bound;thirdfourthq(~ind2,:)];thirdfourthq=thirdfourthq(ind2,:);
firstsecondq=sortrows(firstsecondq,1);thirdfourthq=sortrows(thirdfourthq,1);c=0;i=1;j=1;ind=[];
while(j<=size(firstsecondq,1))%disp(['i-',num2str(i)]);disp(['j-',num2str(j)]);
    if (round(firstsecondq(i,1)*10^(2))/10^(2))==(round(firstsecondq(j,1)*10^(2))/10^(2))&& abs((round(firstsecondq(i,2)*10^(2))/10^(2))-(round(firstsecondq(j,2)*10^(2))/10^(2)))<=3*ste
        if firstsecondq(j,2)>c;c=firstsecondq(j,2);ind=j;
        end;j=j+1;
    else
        new_bound=[new_bound;firstsecondq(ind,:)];i=j;j=i;c=0;
    end
end
c=1000;i=1;j=1;ind=[];
while(j<=size(thirdfourthq,1))%disp(['i-',num2str(i)]);disp(['j-',num2str(j)]);
    if (round(thirdfourthq(i,1)*10^(2))/10^(2))==(round(thirdfourthq(j,1)*10^(2))/10^(2))&& abs((round(thirdfourthq(i,2)*10^(2))/10^(2))-(round(thirdfourthq(j,2)*10^(2))/10^(2)))<=3*ste
        if thirdfourthq(j,2)<c;c=thirdfourthq(j,2);ind=j;
        end
        j=j+1;
    else
        new_bound=[new_bound;thirdfourthq(ind,:)];i=j;j=i;c=1000;
    end
end
%scatter(new_bound(:,1),new_bound(:,2),'filled','c');hold on;

% to compute coordinates and their corresponding angles
store=[];thint=(round(pi/360*10^(sig)))/10^(sig);% 0.5 degrees
for i=1: size(new_bound,1)-1
    x=throat_center_xy(1,1)-new_bound(i,1);y=throat_center_xy(1,2)-new_bound(i,2);r=sqrt(x^2+y^2);
    x_1=throat_center_xy(1,1)-new_bound(i+1,1);y1=throat_center_xy(1,2)-new_bound(i+1,2);r1=sqrt(x_1^2+y1^2);
    the=atan(y/x);
    if (the <=0 && x<0);the= 2*pi+the;end
    if the <0 && y<0 || (the>=0 && x>=0 && y>=0);the= pi+the;end
    the1=atan(y1/x_1);if (the1 <0 && x_1<0);the1= 2*pi+the1;end
    if the1 <0 && y1<0 || (the1>=0 && x_1>=0 && y1>=0);the1= pi+the1;end
    store=[store; new_bound(i,:) r the];
    if i== size(new_bound,1)-1; store=[store; new_bound(end,:) r1 the1];end
end
store=round(store*10^(sig))/10^(sig);
store=unique(store,'rows');%remove duplicate rows
store= sortrows(store,[4,3]);%first sort based on theta and then to break ties sort based on radius
[~, w] = unique( store(:,4), 'stable' );
store(setdiff( 1:numel(store(:,4)), w ),:)=[];%to remove all points with thetas close to each other, those points with higher radii is removed to take into account reentrant angles
%drasticly different r's in both neighbourhood directions and close thetas are outliers
rad=[0;abs(diff(store(:,3)))];%find the peaks of radius difference and remove outliers
ind=find(rad>mean(rad));
store(ind(diff(ind)==1),:)=[];
store=[store [0;diff(store(:,4))]];%add the adjacent difference in theta value as a 4th column
store=[store;store(1,:)];store_new=store;

if vis==1;plot(store_new(:,1),store_new(:,2),'g','Linewidth',3); hold on;
    scatter(throat_center_xy(1,1),throat_center_xy(1,2),'r','filled');hold on;
end

store_new=[store_new linspace(1,size(store_new,1),size(store_new,1)).'];siz= size(store_new,1);
for i=1:siz; store_new(i,4)=(round(store_new(i,4)*10^(sig)))/10^(sig);end% to get all angles in 3 significant digits
for i=1: siz-1
    if i== siz-1
        x=throat_center_xy(1,1)-store_new(i,1);y=throat_center_xy(1,2)-store_new(i,2);r=sqrt(x^2+y^2);
        x_1=throat_center_xy(1,1)-store_new(i+1,1);y1=throat_center_xy(1,2)-store_new(i+1,2);r1=sqrt(x_1^2+y1^2);
        the= store_new(i,4);
        if (round(the*10^(sig)))/10^(sig)>= 2*(round(pi*10^(sig)))/10^(sig); the2=(round(the*10^(sig)))/10^(sig)-2*(round(pi*10^(sig)))/10^(sig);the1=store_new(i+1,4);
        else;the2=the;the1=(round(store_new(i+1,4)*10^(sig)))/10^(sig)+2*(round(pi*10^(sig)))/10^(sig);
        end
        if (round((abs(the2-the1))*10^(sig)))/10^(sig)> (round(thint*10^(sig)))/10^(sig);r_original = [r, r1];angle = [the2 , the1];angle_tointerpolate = the1-mod(abs(the2-the1),thint);
            if angle_tointerpolate ~= the1
                new_r = interp1(angle, r_original, angle_tointerpolate, 'linear','extrap');
                [x2,y2] = pol2cart(angle_tointerpolate,new_r);x2=x2+throat_center_xy(1,1); y2= y2+throat_center_xy(1,2);
                store_new(i+1,1:5)=[x2 y2 new_r angle_tointerpolate 100];
            else
                new_r=r1;
            end
            num_points_to_add= floor((round(((abs(angle_tointerpolate-the2)/thint)-1.00)*10^(sig)))/10^(sig));
            if num_points_to_add~=0
                r_original=[r, new_r];angle=[the2, angle_tointerpolate];a_to_interpolate2=linspace(the2, angle_tointerpolate,num_points_to_add+2);
                new_r = interp1(angle, r_original, a_to_interpolate2(2:end-1), 'linear','extrap');
                [x2,y2] = pol2cart(a_to_interpolate2(2:end-1),new_r);x2=x2+throat_center_xy(1,1); y2= y2+throat_center_xy(1,2);
                store_new=[store_new; x2.' y2.' new_r.' a_to_interpolate2(2:end-1).' repmat(100,1,size(x2,2)).' repmat((i+0.5),1,size(x2,2)).' ];
            end
        end
        break
    end
    x=throat_center_xy(1,1)-store_new(i,1);y=throat_center_xy(1,2)-store_new(i,2);r=sqrt(x^2+y^2);
    x_1=throat_center_xy(1,1)-store_new(i+1,1);y1=throat_center_xy(1,2)-store_new(i+1,2);r1=sqrt(x_1^2+y1^2);
    the= store_new(i,4);the1=store_new(i+1,4);
    if double(abs(the-the1))> double(thint)
        r_original = [r, r1];angle = [the , the1];angle_tointerpolate = the1-mod(abs(the-the1),thint);
        if angle_tointerpolate ~= the1
            new_r = interp1(angle, r_original, angle_tointerpolate, 'linear','extrap');
            [x2,y2] = pol2cart(angle_tointerpolate,new_r);x2=x2+throat_center_xy(1,1); y2= y2+throat_center_xy(1,2);
            store_new(i+1,1:5)=[x2 y2 new_r angle_tointerpolate 100];
        else;new_r=r1;
        end
        num_points_to_add= floor((round(((abs(angle_tointerpolate-the)/thint)-1.00)*10^(sig)))/10^(sig));
        if num_points_to_add~=0
            r_original=[r, new_r];angle=[the, angle_tointerpolate];a_to_interpolate2=linspace(the, angle_tointerpolate,num_points_to_add+2);
            new_r = interp1(angle, r_original, a_to_interpolate2(2:end-1), 'linear','extrap');
            [x2,y2] = pol2cart(a_to_interpolate2(2:end-1),new_r);x2=x2+throat_center_xy(1,1); y2= y2+throat_center_xy(1,2);
            store_new=[store_new; x2.' y2.' new_r.' a_to_interpolate2(2:end-1).' repmat(100,1,size(x2,2)).' repmat((i+0.5),1,size(x2,2)).' ];
        end
    end
end
store_new=sortrows(store_new,4);%isequal(store_new(:,4),(round(linspace(0,store_new(end,4),size(store_new,1))*10^(sig)))/10^(sig))
%idx = setdiff(round(store_new(:,1:2),1),round(surface_profile_xy1,1),'rows');%to find dissimilar parts between boundary and throat pore boundary
if vis ==1
    scatter(store_new(:,1),store_new(:,2),'r', 'filled');hold on;
    title('initial and filtered boundary')
end
%----------------------------------------------------------------------------------END CODE------------------------------------------------------------------------
