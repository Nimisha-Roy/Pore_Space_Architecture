function [store_new,throat_wall_normalized_roughness,throat_center_new]=step5ef_filter(x1,store_new,thint,throat_wall_normalized_roughness,throat_center_new,throat_center_xy1,vis)
% This function filters meso and micro features from input boundary and computes throat wall normalized roughness. It is step 5f of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
radius= store_new(:,3);angle= store_new(:,4);radius(end,1)= radius(1,1);rapprox=[];rapprox_m=[];
ar=[angle radius];ar= sortrows(ar);ar(end,2)= ar(1,2); p1p=ar(1,2); ar(:,2)=ar(:,2)-p1p;radius= ar(:,2);angle= ar(:,1);% finding the radii and angles of equally spaced boundary points
N= size(radius,1);T= thint*(N-1);angle1 = linspace(angle(1,1),T,N);n=N;dx = thint;x =(0:1:n-1).'*dx;X_mags=fft(radius); % number of profile points,spacing in mm,generate x axis data,perform FFT
j=(2:1:floor(n/2)+1); lambda=n*dx./(j-1)';% generate the wavelength array
amp=abs(X_mags(2:floor(n/2)+1,1));har = (1:1:n/2).';har= (1./har);

% cutoff
amp= amp*2/n;rough_cut= 1.2* max(abs(amp))/100;meso_cut= 15* max(abs(amp))/100;amp=amp*n/2;meso_cut_amp= meso_cut*n/2;rough_cut_amp= rough_cut *n/2;
for i= size(har,1):-1:1
    if amp(i,1)> rough_cut_amp;lambdac_r= lambda(i,1);
        break;
    else
        if i ==1;IDX = knnsearch(har,rough_cut_amp); lambdac_r = lambda(IDX,1);
        end
    end
end
for i= size(har,1):-1:1
    if amp(i,1)> meso_cut_amp;lambdac_m= lambda(i,1);
        break;
    else
        if i == 1;IDX = knnsearch(har,meso_cut_amp); lambdac_m = lambda(IDX,1);
        end
    end
end
% removing and measuring roughness of boundary
lambdac= lambdac_r;dx=thint; x2=(-lambdac:dx:lambdac)';alpha=0.4697;S=(1/(alpha*lambdac)).*exp(-pi*(x2/(alpha*lambdac)).^2);%weighing function of gaussian regression filter
S=S/sum(S);% generate the Gaussian filter
m=size(S,1);l=n*thint;n=l/dx;S=[zeros(floor(n/2)-floor(m/2),1); S; zeros(floor(n/2)-floor(m/2)-1,1)]; % length of Gaussian filter; % length of a profile; % number of profile points;% center the filter and zero-pad
Sf=fft(S);j=(2:1:floor(n/2)+1)';wave=n*dx./(j-1); % DFT of S;generate wavelength array for X axis of;transmission plot

% low pass filter specified wavelength using gaussian regression filter of second order for 50% amplitude attenuation at cutoff
const = sqrt(log(2)/2/pi/pi);
for k = 1:n;p = (1:1:n)'; % For each position k (center of the filter), generate the filter over the entire profile length
    S = (1/sqrt(2*pi*pi)/const/lambdac).*exp(-0.5*((k-p)*dx/const/lambdac).^2);% generate weighting function
    S = S/sum(S); % normalize to unit sum  and three variables given in Eq. 9.9. Ignore the dx term in Eq. 9.9 because it is a constant.
    x= (k-p)*dx;M=[];Q=[];M(1,1) = sum(S.*x.^4);M(1,2) = sum(S.*x.^3);M(1,3) = sum(S.*x.^2);M(2,1) = M(1,2);M(2,2) = M(1,3);M(2,3) = sum(S.*x);
    M(3,1) = M(1,3);M(3,2) = M(2,3);M(3,3) = sum(S);
    Q(1,1) = sum(radius.*S.*x.^2);Q(2,1) = sum(radius.*S.*x);Q(3,1) = sum(radius.*S);
    P = inv(M) * Q; % determine array P containing the values of % three unknown
    A = P(1); B = P(2); C = P(3);
    rapprox(k,1) = C; % determine roughness free profile value at that location
end

% reconstruction
radius= radius + p1p;rapprox= rapprox + p1p;angle1=angle1.';x_new=[];y_new=[];
for i=1:size(angle1,1)
    x_new=[x_new; rapprox(i,1)*cos(angle1(i,1))];y_new=[y_new; rapprox(i,1)*sin(angle1(i,1))];
end
x_new= x_new+ throat_center_xy1(1,1);y_new= y_new+ throat_center_xy1(1,2);x_new(end,1)=x_new(1,1);y_new(end,1)=y_new(1,1);x_old=[];y_old=[];
for i=1:size(angle1,1)
    x_old=[x_old; radius(i,1)*cos(angle1(i,1))];y_old=[y_old; radius(i,1)*sin(angle1(i,1))];
end
x_old= x_old+ throat_center_xy1(1,1);y_old= y_old+ throat_center_xy1(1,2);
mean_dif= mean(abs(rapprox - radius));
throat_wall_normalized_roughness=[throat_wall_normalized_roughness;sqrt(mean((abs(rapprox - radius)- mean_dif).^2))];
throat_center_new(x1,7)=sqrt(mean((abs(rapprox - radius)- mean_dif).^2));

% FOR REMOVING MESO-SCALE FEATURES
radius= radius - p1p;rapprox=rapprox-p1p;lambdac= lambdac_m;dx=thint;x2=(-lambdac:dx:lambdac)';alpha=0.4697;S=(1/(alpha*lambdac)).*exp(-pi*(x2/(alpha*lambdac)).^2);
S=S/sum(S);m=size(S,1);l=n*thint; n=l/dx;S=[zeros(floor(n/2)-floor(m/2),1); S; zeros(floor(n/2)-floor(m/2)-1,1)];
Sf=fft(S); j=(2:1:floor(n/2)+1)';wave=n*dx./(j-1); const = sqrt(log(2)/2/pi/pi);

for k = 1:n;p =(1:1:n)'; % For each position k (center of the filter), generate the filter over the entire profile length
    S = (1/sqrt(2*pi*pi)/const/lambdac).*exp(-0.5*((k-p)*dx/const/lambdac).^2);S = S/sum(S);x= (k-p)*dx;M=[];Q=[];
    M(1,1) = sum(S.*x.^4);M(1,2) = sum(S.*x.^3);M(1,3) = sum(S.*x.^2);M(2,1) = M(1,2);M(2,2) = M(1,3);M(2,3) = sum(S.*x);
    M(3,1) = M(1,3);M(3,2) = M(2,3);M(3,3) = sum(S);
    Q(1,1) = sum(radius.*S.*x.^2);Q(2,1) = sum(radius.*S.*x);Q(3,1) = sum(radius.*S);
    P = inv(M) * Q;A = P(1); B = P(2); C = P(3);
    rapprox_m(k,1) = C;
end
%figure, plot(xxx,radius,'k','LineWidth', 1);hold on;plot(xxx,rapprox,'r','LineWidth', 1.5);hold on;
%plot(x,radius,'k','LineWidth', 1);hold on;plot(x,rapprox_m,'g','LineWidth', 1.5);xlabel('Angle (radians)');ylabel('Radial distance (pixels)');title('Original,roughness removed and meso features removed profiles');

% reconstruction
radius= radius + p1p;rapprox_m= rapprox_m + p1p;x_new1=[];y_new1=[];
for i=1:size(angle1,1)
    x_new1=[x_new1; rapprox_m(i,1)*cos(angle1(i,1))];y_new1=[y_new1; rapprox_m(i,1)*sin(angle1(i,1))];
end
x_new1=  x_new1+ throat_center_xy1(1,1);y_new1= y_new1+ throat_center_xy1(1,2);x_new1(end,1)=x_new1(1,1);y_new1(end,1)=y_new1(1,1);x_old1=[];y_old1=[];
for i=1:size(angle1,1)
    x_old1=[x_old1; radius(i,1)*cos(angle1(i,1))];y_old1=[y_old1; radius(i,1)*sin(angle1(i,1))];
end
x_old1= x_old1+ throat_center_xy1(1,1);y_old1= y_old1+ throat_center_xy1(1,2);
if vis==1
    figure,plot(x_old(:,1), y_old(:,1), 'k', 'LineWidth', 2);hold on
    plot(x_new1(:,1), y_new1(:,1), 'g', 'LineWidth', 2);title('original,roughness removed and meso feature removed boundary plot');
end
store_new(:,1:2)=[x_new1 y_new1];
%--------------------------------------------------------------------------END CODE---------------------------------------------------------------------------------------