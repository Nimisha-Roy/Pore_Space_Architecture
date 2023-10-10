function locals = step1_local_maxima(vox3d,ste,mean_input_radius,porosity)
% This function computes the center of pore volumes in the input structure. It is part of Step 1 of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
% vox3d is the voxelized structure
vox3d1 = bwdist(~vox3d); % to find distance of each voxel from the nearest zero voxel which is the solid voxel
sb = size(vox3d1);

%Local maxima stability criteria for optimal search block size
if (((1.5*mean_input_radius*porosity)/ste)>1)
    ii=round((1.5*mean_input_radius*porosity)/ste);
else
    ii=5;
end

% We are trying to find the peaks among void pixels using euclidean distance transform. For this purpose, first for each void pixel, we compare its EDT value with its nearest neighbours, which is 24 for a voxel.
% We then keep on increasing the extension of the neighbours to 2nd degree, 3rd degree and so on until the number of obtained peaks stabilizes and we can say that the right peaks are obtained.

nn = 1:50; % nn is the extension rate. nn equals 1 is considered as immediate neighbors and n = 2 is considerd for secondary neighbours too

% The size of the matrix to be used for comparison is important for accuracy of the search. For example, for n = 2 we have a matrix with size 5x5x5,
% in which 124 members of this matrix are the neighbours and one is the voxel which is going to be compared.

step = nn(ii)*2 + 1;
bound = (step -1)/2;% number of members to be compared on each side of the current voxel
main = zeros(sb(1)+(bound*2),sb(2)+(bound*2),sb(3)+(bound*2)); % the size of the main matrix is increased/fabricated by bound size on either sides of each dimension to take into account boundary effects % and it has been assigned solid pixel values to not disturb our comparision while finding peaks.
sm = size(main);
main(nn(ii)+1:end-nn(ii),nn(ii)+1:end-nn(ii),nn(ii)+1:end-nn(ii)) = vox3d1; % assigning EDT values to all voxels of the actual structure in the middle of the fabricated matrix
counter = 1; % the local centers counter
positions = [];
locals = [];

for i = nn(ii)+1:sm(1)-nn(ii) % computation starts from the actual voids, not the fabricated voxels
    for j = nn(ii)+1:sm(2) - nn(ii)
        for k = nn(ii)+1:sm(3) - nn(ii)
            if main(i,j,k) == 0 % to compute and compare only if voids, for solid voxels, b value is 0 and for fabricated voxels, value is zero too
                continue;
            end
            temp = main(i-(bound+1)+1:i-(bound+1)+step,j-(bound+1)+1:j-(bound+1)+step,k-(bound+1)+1:k-(bound+1)+step); % the temp ECD matrix for a void voxel with designated neighbours according to nn
            nonzerosTemp = nonzeros(temp); MAX = max(nonzerosTemp); MIN = 0;%min(min(min(temp)));
            NumberFZ = sum(nonzerosTemp-MAX == 0); % to check the number of max values , if it is greater than 1
            if MAX == temp(nn(ii)+1,nn(ii)+1,nn(ii)+1) && NumberFZ == 1 % if current voxel is the only peak
                positions(counter,:) = [i-bound j-bound k-bound]; % positions of actual structure matrix, vox3d stored
                counter = counter+1;
            elseif MAX == temp(nn(ii)+1,nn(ii)+1,nn(ii)+1) && NumberFZ ~= 1
                main(i,j,k) = MIN;
            end
        end
    end
end

for t = 1:length(positions)
    locals(t,1) = sub2ind(sb,positions(t,1),positions(t,2),positions(t,3) );
end
locals = (sort(locals));
%----------------------------------------------------------------- END CODE--------------------------------------------------------------------------