function COOR_DP=step5g_inflexion_dominant(xy_r,vis)
%This function finds the dominant points along a boundary. It is a sub-function used in inflexion_combined.m and is part of step 5g,h of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
debug=0;
xy_r=[xy_r(:,2),xy_r(:,1)];xy_r=[xy_r;xy_r(1,:)];
x = xy_r(:,2); y = xy_r(:,1); % x and y components of coordinate
g= pdist2(xy_r,xy_r,'euclidean');
nm = max(g(:));
[ro,~] = find(g== nm);
index = ro(1,1);index1 = ro(2,1);
c= sqrt((xy_r(ro(1,1),2)-xy_r(ro(2,1),2))^2 + (xy_r(ro(1,1),1)-xy_r(ro(2,1),1))^2);
n=c*0.01;
if debug ==1
    figure,scatter(xy_r(:,2),xy_r(:,1),'r');h=impoint(gca, xy_r(ro(1,1),2),xy_r(ro(1,1),1));setColor(h,'g');
    h1=impoint(gca, xy_r(ro(2,1),2),xy_r(ro(2,1),1));setColor(h1,'g');line([xy_r(ro(1,1),2) xy_r(ro(2,1),2)],[xy_r(ro(1,1),1) xy_r(ro(2,1),1)], 'Color','b','LineWidth',2);hold on
end
% the boundary is divided and operated on in 2 parts
% for first half
if max(index,index1) <= (length(x))/2 % index is location of the two extreme points in xy_r matrix
    v=[ (xy_r(min(index,index1),2)),(xy_r(min(index,index1),1)),0; (xy_r(max(index,index1),2)),(xy_r(max(index,index1),1)),0 ]; % v stores coordinates of the two extreme points
    im=[ min(index,index1); max(index,index1)]; % index of extreme points are stored in im
    dsp=c; % c is Maximum distance
    [im_new,v_new]=step5g_farthest_points1(dsp,n,im,x,y,v,c,debug);
else
    pp =[];
    ss= max(index,index1)- min(index,index1) +1;
    ssd= length(x)- max(index,index1);
    pp(1:ss,2)= x( min(index,index1): max(index,index1));
    pp(1:ss,1)= y( min(index,index1): max(index,index1));
    pp(ss+1:ss+ssd,2)= x( max(index,index1)+1: length(x));
    pp(ss+1:ss+ssd,1)= y( max(index,index1)+1: length(x));
    pp(ss+ssd+1:length(x),2)= x( 1: min(index,index1)-1);
    pp(ss+ssd+1:length(x),1)= y( 1: min(index,index1)-1);
    v=[ (pp(1,2)),(pp(1,1)),0; (pp(ss,2)),(pp(ss,1)),0 ];
    im=[ 1; ss];
    dsp=c;
    [im_new,v_new]=step5g_farthest_points2(dsp,n,v,pp,c,im,debug);
end
% for 2nd half
if max(index,index1) <= (length(x))/2
    vn=[ (xy_r(max(index,index1),2)),(xy_r(max(index,index1),1)),0; (xy_r(min(index,index1),2)),(xy_r(min(index,index1),1)),0 ];
    imn=[ max(index,index1); length(x)];
    dspn=c;
    [im_newn,v_newn]=step5g_farthest_points1(dspn,n,imn,x,y,vn,c,debug);
else
    pp =[];
    ss= max(index,index1)- min(index,index1) +1;
    ssd= length(x)- max(index,index1);
    pp(1:ss,2)= x( min(index,index1): max(index,index1));
    pp(1:ss,1)= y( min(index,index1): max(index,index1));
    pp(ss+1:ss+ssd,2)= x( max(index,index1)+1: length(x));
    pp(ss+1:ss+ssd,1)= y( max(index,index1)+1: length(x));
    pp(ss+ssd+1:length(x),2)= x( 1: min(index,index1)-1);
    pp(ss+ssd+1:length(x),1)= y( 1: min(index,index1)-1);
    vn=[ (pp(ss,2)),(pp(ss,1)),0; (pp(1,2)),(pp(1,1)),0 ];
    imn=[ ss+1; length(pp)];
    dspn=c;
    [im_newn,v_newn]=step5g_farthest_points2(dspn,n,vn,pp,c,imn,debug);
end
v_new(:,3) = 0;v_newn(:,3) = 0;
INDEX1= [im_new; im_newn];
COOR_DP1=[v_new ; v_newn];
COOR_DP=[];INDEX= [];
po=1;
for i= 1:size(COOR_DP1,1)
    if i== size(COOR_DP1,1)
        COOR_DP(po,1:2)= COOR_DP1(i,1:2);
        INDEX(po,:)= INDEX1(i,:);
        break
    end
    if COOR_DP1(i,1) ~= COOR_DP1(i+1,1) || COOR_DP1(i,2) ~= COOR_DP1(i+1,2)
        COOR_DP(po,1:2)= COOR_DP1(i,1:2);
        INDEX(po,:)= INDEX1(i,:);
        po= po+1;
    end
end
if vis==1
    scatter(COOR_DP(:,1),COOR_DP(:,2),'m','filled')
    hold on
end
%------------------------------------------------------------------------------END CODE------------------------------------------------------------------------------
