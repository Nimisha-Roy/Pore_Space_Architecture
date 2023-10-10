function [store_new_combined]= step5ef_filtered_combined(combined_profile,ste,sig,thint,vis)
%This function filters the combined boundary of projected throat and adjacent body voxels. It is part of Step 5e f of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
combined_profile=(round(combined_profile*10^(sig)))/10^(sig);combined_profile=unique(combined_profile,'rows');
bound = boundary(combined_profile(:,1),combined_profile(:,2),1);% approximate boundary point coordinates to emsure boundary coordinates plus a bunch of infiltrated coordinates
bound= unique(bound);new_bound=[combined_profile(bound,1),combined_profile(bound,2)];
throat_center_xy= [mean(combined_profile(:,1)),mean(combined_profile(:,2))];

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
firstsecondq=sortrows(firstsecondq,1);thirdfourthq=sortrows(thirdfourthq,1);c=0;i=1;j=1;ind=j;

while(j<=size(firstsecondq,1))%disp(['i-',num2str(i)]);disp(['j-',num2str(j)]);
    if (round(firstsecondq(i,1)*10^(2))/10^(2))==(round(firstsecondq(j,1)*10^(2))/10^(2))&& abs((round(firstsecondq(i,2)*10^(2))/10^(2))-(round(firstsecondq(j,2)*10^(2))/10^(2)))<=3*ste
        if firstsecondq(j,2)>c;c=firstsecondq(j,2);ind=j;
        end
        j=j+1;
    else;new_bound=[new_bound;firstsecondq(ind,:)];i=j;j=i;c=0;
    end
end
c=1000;i=1;j=1;ind=j;
while(j<=size(thirdfourthq,1))%disp(['i-',num2str(i)]);disp(['j-',num2str(j)]);
    if (round(thirdfourthq(i,1)*10^(2))/10^(2))==(round(thirdfourthq(j,1)*10^(2))/10^(2))&& abs((round(thirdfourthq(i,2)*10^(2))/10^(2))-(round(thirdfourthq(j,2)*10^(2))/10^(2)))<=3*ste
        if thirdfourthq(j,2)<c;c=thirdfourthq(j,2);ind=j;
        end
        j=j+1;
    else
        new_bound=[new_bound;thirdfourthq(ind,:)];i=j;j=i;c=1000;
    end
end
if vis==1
    figure,scatter(new_bound(:,1),new_bound(:,2),'c');hold on;
end
% to compute coordinates and their corresponding angles
store=[];
for i=1: size(new_bound,1)-1
    x=throat_center_xy(1,1)-new_bound(i,1);y=throat_center_xy(1,2)-new_bound(i,2);r=sqrt(x^2+y^2);
    x_1=throat_center_xy(1,1)-new_bound(i+1,1);y1=throat_center_xy(1,2)-new_bound(i+1,2);r1=sqrt(x_1^2+y1^2);
    the=atan(y/x);
    if (the <=0 && x<0)
        the= 2*pi+the;
    end
    if the <0 && y<0 || (the>=0 && x>=0 && y>=0)
        the= pi+the;
    end
    the1=atan(y1/x_1);
    if (the1 <0 && x_1<0)
        the1= 2*pi+the1;
    end
    if the1 <0 && y1<0 || (the1>=0 && x_1>=0 && y1>=0)
        the1= pi+the1;
    end
    store=[store; new_bound(i,:) r the];
    if i== size(new_bound,1)-1
        store=[store; new_bound(end,:) r1 the1];
    end
end
store=round(store*10^(sig))/10^(sig);
store=unique(store,'rows');%remove duplicate rows
store= sortrows(store,[4,3]);%first sort based on theta and then to break ties sort based on radius
[~, w] = unique( store(:,4), 'stable' );
store(setdiff( 1:numel(store(:,4)), w ),:)=[];%to remove all points with thetas close to each other, those points with higher radii is removed to take into account reentrant angles
%drasticly different r's in both neighbourhood directions and close thetas are outliers,remove them
rad=[0;abs(diff(store(:,3)))];%find the peaks of radius difference and remove outliers
ind=find(rad>mean(rad));
store(ind(diff(ind)==1),:)=[];
store=[store [0;diff(store(:,4))]];%add the adjacent difference in theta value as a 4th column
store=[store;store(1,:)];store_new_combined=store;

if vis==1
    plot(store_new_combined(:,1),store_new_combined(:,2),'g','Linewidth',3); hold on;
    scatter(throat_center_xy(1,1),throat_center_xy(1,2),'r','filled');hold on;
end

%add points to intervals higher than 0.5 degrees
store_new_combined=[store_new_combined linspace(1,size(store_new_combined,1),size(store_new_combined,1)).'];siz= size(store_new_combined,1);
for i=1: siz-1
    if i== siz-1
        x=throat_center_xy(1,1)-store_new_combined(i,1);y=throat_center_xy(1,2)-store_new_combined(i,2);r=sqrt(x^2+y^2);
        x_1=throat_center_xy(1,1)-store_new_combined(i+1,1);y1=throat_center_xy(1,2)-store_new_combined(i+1,2);r1=sqrt(x_1^2+y1^2);
        the= store_new_combined(i,4);
        if (round(the*10^(sig)))/10^(sig)>= 2*(round(pi*10^(sig)))/10^(sig)
            the2=(round(the*10^(sig)))/10^(sig)-2*(round(pi*10^(sig)))/10^(sig);the1=store_new_combined(i+1,4);
        else
            the2=the;the1=(round(store_new_combined(i+1,4)*10^(sig)))/10^(sig)+2*(round(pi*10^(sig)))/10^(sig);
        end
        if (round((abs(the2-the1))*10^(sig)))/10^(sig)> (round(thint*10^(sig)))/10^(sig);r_original = [r, r1];angle = [the2 , the1];angle_tointerpolate = the1-mod(abs(the2-the1),thint);
            if angle_tointerpolate ~= the1
                new_r = interp1(angle, r_original, angle_tointerpolate, 'linear','extrap');
                [x2,y2] = pol2cart(angle_tointerpolate,new_r);x2=x2+throat_center_xy(1,1); y2= y2+throat_center_xy(1,2);
                store_new_combined(i+1,1:5)=[x2 y2 new_r angle_tointerpolate 100];
            else
                new_r=r1;
            end
            num_points_to_add= floor((round(((abs(angle_tointerpolate-the2)/thint)-1.00)*10^(sig)))/10^(sig));
            if num_points_to_add~=0
                r_original=[r, new_r];angle=[the2, angle_tointerpolate];a_to_interpolate2=linspace(the2, angle_tointerpolate,num_points_to_add+2);
                new_r = interp1(angle, r_original, a_to_interpolate2(2:end-1), 'linear','extrap');
                [x2,y2] = pol2cart(a_to_interpolate2(2:end-1),new_r);x2=x2+throat_center_xy(1,1); y2= y2+throat_center_xy(1,2);
                store_new_combined=[store_new_combined; x2.' y2.' new_r.' a_to_interpolate2(2:end-1).' repmat(100,1,size(x2,2)).' repmat((i+0.5),1,size(x2,2)).' ];
            end
        end
        break;
    end
    x=throat_center_xy(1,1)-store_new_combined(i,1);y=throat_center_xy(1,2)-store_new_combined(i,2);r=sqrt(x^2+y^2);
    x_1=throat_center_xy(1,1)-store_new_combined(i+1,1);y1=throat_center_xy(1,2)-store_new_combined(i+1,2);r1=sqrt(x_1^2+y1^2);
    the= store_new_combined(i,4);the1=store_new_combined(i+1,4);
    if double(abs(the-the1))> double(thint);r_original = [r, r1];angle = [the , the1];angle_tointerpolate = the1-mod(abs(the-the1),thint);
        if angle_tointerpolate ~= the1
            new_r = interp1(angle, r_original, angle_tointerpolate, 'linear','extrap');
            [x2,y2] = pol2cart(angle_tointerpolate,new_r);x2=x2+throat_center_xy(1,1); y2= y2+throat_center_xy(1,2);
            store_new_combined(i+1,1:5)=[x2 y2 new_r angle_tointerpolate 100];
        else;new_r=r1;
        end;num_points_to_add= floor((round(((abs(angle_tointerpolate-the)/thint)-1.00)*10^(sig)))/10^(sig));
        if num_points_to_add~=0
            r_original=[r, new_r];angle=[the, angle_tointerpolate];a_to_interpolate2=linspace(the, angle_tointerpolate,num_points_to_add+2);
            new_r = interp1(angle, r_original, a_to_interpolate2(2:end-1), 'linear','extrap');
            [x2,y2] = pol2cart(a_to_interpolate2(2:end-1),new_r);x2=x2+throat_center_xy(1,1); y2= y2+throat_center_xy(1,2);
            store_new_combined=[store_new_combined; x2.' y2.' new_r.' a_to_interpolate2(2:end-1).' repmat(100,1,size(x2,2)).' repmat((i+0.5),1,size(x2,2)).' ];
        end
    end
end
store_new_combined=sortrows(store_new_combined,4);
%isequal(store_new_combined(:,4),(round(linspace(0,store_new_combined(end,4),size(store_new_combined,1))*10^(sig)))/10^(sig))
%figure, scatter(mid_profile_xy1(:,1),mid_profile_xy1(:,2),'b');hold on;scatter(body_center_xy1(:,1),body_center_xy1(:,2),'r','filled');hold on;
%scatter(body1_profile_xy1(:,1),body1_profile_xy1(:,2),'g');hold on;scatter(body2_profile_xy1(:,1),body2_profile_xy1(:,2),'g');hold on;

% filtering meso features from combined boundary
radius= store_new_combined(:,3);angle= store_new_combined(:,4);radius(end,1)= radius(1,1);rapprox_m=[];
ar=[angle radius];ar= sortrows(ar);ar(end,2)= ar(1,2); p1p=ar(1,2); ar(:,2)=ar(:,2)-p1p;radius= ar(:,2);angle= ar(:,1);% finding the radii and angles of equally spaced boundary points
N= size(radius,1);T= thint*(N-1);angle1 = linspace(angle(1,1),T,N);n=N;dx = thint;x =(0:1:n-1).'*dx;X_mags=fft(radius); % number of profile points,spacing in mm,generate x axis data,perform FFT
j=(2:1:floor(n/2)+1); lambda=n*dx./(j-1)';% generate the wavelength array
amp=abs(X_mags(2:floor(n/2)+1,1));har = (1:1:n/2).';har= (1./har);
% cutoff
amp= amp*2/n;meso_cut= 15* max(abs(amp))/100;amp=amp*n/2;meso_cut_amp= meso_cut*n/2;
for i= size(har,1):-1:1
    if amp(i,1)> meso_cut_amp
        harc_m= har(i,1);lambdac_m= lambda(i,1);
        break;
    else
        if i == 1
            IDX = knnsearch(har,meso_cut_amp);
            harc_m= har(IDX,1);lambdac_m = lambda(IDX,1);
        end
    end
end

lambdac= lambdac_m;dx=thint;x2=(-lambdac:dx:lambdac)';alpha=0.4697;S=(1/(alpha*lambdac)).*exp(-pi*(x2/(alpha*lambdac)).^2);
S=S/sum(S);m=size(S,1);l=n*thint; n=l/dx;S=[zeros(floor(n/2)-floor(m/2),1); S; zeros(floor(n/2)-floor(m/2)-1,1)];
Sf=fft(S); j=(2:1:floor(n/2)+1)';wave=n*dx./(j-1); const = sqrt(log(2)/2/pi/pi);

for k = 1:n;p =(1:1:n)'; % For each position k (center of the filter), generate the filter over the entire profile length
    S = (1/sqrt(2*pi*pi)/const/lambdac).*exp(-0.5*((k-p)*dx/const/lambdac).^2);S = S/sum(S);x= (k-p)*dx;M=[];Q=[];
    M(1,1) = sum(S.*x.^4);M(1,2) = sum(S.*x.^3);M(1,3) = sum(S.*x.^2);M(2,1) = M(1,2);M(2,2) = M(1,3);M(2,3) = sum(S.*x);
    M(3,1) = M(1,3);M(3,2) = M(2,3);M(3,3) = sum(S);
    Q(1,1) = sum(radius.*S.*x.^2);Q(2,1) = sum(radius.*S.*x);Q(3,1) = sum(radius.*S);
    P = inv(M) * Q;A = P(1); B = P(2); C = P(3);rapprox_m(k,1) = C;
end

% reconstruction
radius= radius + p1p;rapprox_m= rapprox_m + p1p;x_new1=[];y_new1=[];angle1=angle1.';
for i=1:size(angle1,1)
    x_new1=[x_new1; rapprox_m(i,1)*cos(angle1(i,1))];y_new1=[y_new1; rapprox_m(i,1)*sin(angle1(i,1))];
end
x_new1=  x_new1+ throat_center_xy(1,1);y_new1= y_new1+ throat_center_xy(1,2);x_new1(end,1)=x_new1(1,1);y_new1(end,1)=y_new1(1,1);x_old1=[];y_old1=[];
for i=1:size(angle1,1)
    x_old1=[x_old1; radius(i,1)*cos(angle1(i,1))];y_old1=[y_old1; radius(i,1)*sin(angle1(i,1))];
end
x_old1= x_old1+ throat_center_xy(1,1);y_old1= y_old1+ throat_center_xy(1,2);
store_new_combined(:,1:2)=[x_new1 y_new1];
if vis==1
    plot(store_new_combined(:,1),store_new_combined(:,2),'r','LineWidth',2);hold on;
    title('Unfiltered and filtered combined boundary')
end
%---------------------------------------------------------------------END CODE----------------------------------------------------------------------------------------
