%% This software segments the pore phase of a binary input microstructure into an interconnected network of bulky pore bodies and narrow pore throats, and computes various individual and collective properties of the pore space.%%
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Give an input binary 3d structure (solid phase voxels are designated as 0 and pore phase voxels are designated as 1), the code outputs the following properties of the structure
% 3D structure segmented into particle phase, individual pore bodies and individual pore throats
% Sizes and shapes of individual pore bodies
% Sizes and shapes of individual pore throats
% Sizes of individual constrictions
% Sharpness along the lengths of individual pore throats
% Wall roughness of individual pore throats
% Geometrical tortuosity of the structure in a given direction
% Coordination number of pore bodies and pore throats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INPUTS (EXAMPLE INPUTS HERE- EDIT THIS) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname ='test_input.mat';%name of input matfile- a 3d structure comprising solid phase and void phase labelled as 0 and 1 respectively 
res = 1/40;% resolution of the image => in voxels per mm. A resolution of 1/40 means 40 voxels are present in 1 mm length of the image. 
% A recommended resolution is include about 40 voxels along the length of 1 particle for error free computations
mpd = 1;%mean particle diameter in mm
output_fname ='test_output.mat';% output matfilename to store results
tortuous_axis = 3;
vtk_file = 'labelled_test.vtk';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MAIN ALGORITHM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_pore_body, num_throat, labelled_output, pore_body_sizes, pore_throat_sizes, constriction_sizes, lsi, throat_wall_normalized_roughness, tortuosity]=pore_space_architecture(fname, res, output_fname, mpd, tortuous_axis, vtk_file);

% Parameters:

% fname (string) – filename of Input binary structure 
% res (int) – resolution of the input structure
% output_fname (string) – output mat filename to store all results
% mpd (int) – mean particle diameter of input structure
% tortuous_axis (int, optional) - direction along which tortuosity needs to be computed, if tortuosity is needed as an output. 1 - along x axis, 2 - along y axis, 3 - along z axis. Default - 1
% vtk_file (string, optional) - filename to save vtk output of labelled structure if needed. Default - False

% Returns:

% num_pore_body (int) - number of pore bodies in the input structure
% num_throat (int) - number of pore throats in the input structure
% labelled_output (3D array) – labelled input structure, 0- particle phase, 1 to num_pore_body - Individual pore bodies, num_pore_body + 1 to num_pore_body + num_throat - Individual pore throats
% pore body sizes (2D array) - lengths of axes a,b,c of the ellipitic fit to the individual pore bodies
% pore throat sizes (2D array) - lengths of axes a,b,c of the ellipitic fit to the individual pore throats
% constriction sizes (1D array) - lengths of individual constriction surfaces
% lsi (1D array) - longitudinal sharpness indices of individual throats
% throat_wall_normalized_roughness (1D array) - normalized roughness of individual throat walls in percentage
% tortuosity (1D array) - tortuosity of paths along specified direction