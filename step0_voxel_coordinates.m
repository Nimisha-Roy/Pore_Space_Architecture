function voxels=step0_voxel_coordinates(X_start,Y_start,X_end,Y_end,Z_start,Z_end,res)
% This function returns the coordinates of each voxel of the input structure. This is a preprocessing step of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
[X_grid,Y_grid,Z_grid] = ndgrid(X_start:res:X_end,Y_start:res:Y_end,Z_start:res:Z_end);
lTI = size(X_grid,1);lTJ = size(Y_grid,2);lTK = size(Z_grid,3);

% Finding centroid of the mesh because ndgrid gives the coordinates of the vertices of the grid and should be changed to the centroid of each cubic voxel
TI = X_grid(1:lTI-1,1:lTJ-1,1:lTK-1) + res/2 ; sTI = size(TI); TI = (reshape(TI,numel(TI),1));clearvars X_grid;
TJ = Y_grid(1:lTI-1,1:lTJ-1,1:lTK-1) + res/2 ; TJ = (reshape(TJ,numel(TJ),1));clearvars Y_grid;
TK = Z_grid(1:lTI-1,1:lTJ-1,1:lTK-1) + res/2 ; TK = (reshape(TK,numel(TK),1));clearvars Z_grid;
voxels = [TI,TJ,TK]; clearvars TI TJ TK
%----------------------------------------------------------------- END CODE--------------------------------------------------------------------------
