function [num_pore_body, num_throat, labelled_output, pore_body_sizes, pore_throat_sizes, constriction_length, lsi, throat_wall_normalized_roughness,tortuosity]=pore_space_architecture(fname, res, output_fname, mpd, varargin)
%This is the main function ensuing all five steps of the algorithm
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
% Parameters:

% fname (string) – filename of Input binary structure
% res (int) – resolution of the input structure
% output_fname (string) – output mat filename to store all results
% mpd (int) – mean particle diameter of input structure in mm
% tortuous_axis (int, optional) - direction along which tortuosity needs to be computed, if tortuosity is needed as an output. 1 - along x axis, 2 - along y axis, 3 - along z axis. Default - 1
% vtk_file (string, optional) - filename to save vtk output of labelled structure if needed. Default - False

% Returns:

% num_pore_body (int) - number of pore bodies in the input structure
% num_throat (int) - number of pore throats in the input structure
% labelled_output (3D array) – labelled input structure, 0- particle phase, 1 to num_pore_body - Individual pore bodies, num_pore_body + 1 to num_pore_body + num_throat - Individual pore throats
% pore body sizes (2D array) - lengths of axes a,b,c of the ellipitic fit to the individual pore bodies in mm
% pore throat sizes (2D array) - lengths of axes a,b,c of the ellipitic fit to the individual pore throats in mm
% constriction length (1D array) - lengths of individual constriction surfaces in mm
% lsi (1D array) - longitudinal sharpness indices of individual throats
% throat_wall_normalized_roughness (1D array) - normalized roughness of individual throat walls in percentage
% tortuosity (1D array) - tortuosity of paths along specified direction

%----------------------------------------------------------------- START CODE--------------------------------------------------------------------------

%************************************************************* PREPROCESSING STEPS %********************************************************************

%Read input structure
load(fname,'vox3d');%loading the input structure

%get the size and limits of input structure in cartesian coordinates
size_input = size(vox3d);
X_start = 0;  Y_start = 0; Z_start = 0;
X_end = size_input(1)*res; Y_end = size_input(2)*res; Z_end = size_input(3)*res;

%computing porosity of the input structure
porosity = sum(vox3d(:))/(size_input(1)*size_input(2)*size_input(3));

% function to obtain the cartesian coordinates of each voxel
voxels = step0_voxel_coordinates(X_start,Y_start,X_end,Y_end,Z_start,Z_end,res);

%************************************************ STEP 1: PORE VOLUMES BASED ON EDT AND MEDA %************************************************************

h = waitbar(0,"Starting computation");
% function to compute the local maxima of pores to compute centers of pores
locals = step1_local_maxima(vox3d,res, mpd/2, porosity);
waitbar(0.15,h,"Local Maxima computed");

% function to separate pores into pore volumes
[separator,farthest_voxel] = step1_separate_porevolumes(locals,voxels,vox3d);

%******************************* STEPS 2 & 3: CONSTRICTION SURFACES AND PORE THROAT CENTERS BETWEEN ADJACENT PORE VOLUMES %********************************

% function to find neighbours of each pore volume
[neighbour]= step2_neighbours(separator,farthest_voxel,locals,voxels);

num_pore_body= length(locals);
waitbar(0.25,h,"Pore Volumes computed");

% function to find pore throat centers and constriction surfaces between pores
% note that only those volume elements are throats that have constriction size less than min of distance from respective body centers
[throat_center,neighbour,constriction_length]= step3_throatcenters(separator,voxels,res,locals,num_pore_body,neighbour);
waitbar(0.45,h,"Throat centers computed");

%********************************************** STEP 4: PORE BODIES AND THROATS BASED ON MEDA %*************************************************************

% function for segmentation based on throat centers and pore body centers
[throat_body,throat_center,neighbour,locals]= step4_throat_body(locals,voxels,throat_center,vox3d,neighbour,num_pore_body);

% to extract the throat-body surfaces
[throat_body_surface]=step4_throat_body_surface(throat_body,voxels,neighbour);
waitbar(0.60,h,"Throat and body segmentation done");

%****************************** STEP 5: GEOMETRICAL MODIFICATION OF BODY AND THROAT LENGTHS BASED ON BOUNDARY CURVATURES %************************************

% finding actual pore throat and pore body volumes, throat wall roughness and longitudinal sharpness index of each throat based on inflexion points (which is based on shape
% features of particles making up the throat- by removal of micro and meso level features based on fourier signatures and computing inflexion point
% of the macro scale feature along the direction of the line joining the throat centers)

sig=2;%number of significant digits after decimal point based on resolution
sig1=2;%number of significant digits after decimal point based on precision of angle interval to be considered for extracting constant interval boundary of projected throat voxels
num_pore_body= length(locals);
[labelled_output,lsi,throat_wall_normalized_roughness,throat_center_new,locals_new]= step5_main(throat_body,voxels,num_pore_body,throat_center,res,locals,throat_body_surface,sig,sig1,0);% the last parameter defines no visualization
num_pore_body= length(nonzeros(locals_new));
waitbar(0.90,h,"Final segmentation based on modified lengths done");

%************************************************************ POST-PROCESSING %*********************************************************************************

% compute pore body and throat sizes
indices = unique(nonzeros(labelled_output));
num_throat = size(indices,1)-num_pore_body;
[size_body,size_throat]=step6_feature_lengths(labelled_output, num_pore_body,voxels,X_start,X_end, Y_start,Y_end, Z_start,Z_end, indices);
pore_body_sizes = size_body(~all(roundn(size_body,-3)==0,2),:);
pore_throat_sizes = size_throat(~all(roundn(size_throat,-3)==0,2),:);

%compute tortuosity
if varargin{1}
    tortuous_axis= varargin{1};
    [~, unconnected_pores]=step6_connectivity(neighbour, labelled_output);
    [tortuosity,~]=step6_tortuosity(labelled_output,throat_center_new,unconnected_pores,voxels,num_pore_body,tortuous_axis,locals_new,size_input,neighbour,res);
end

%save vtk files of the segmented structure
if varargin {2}
    vtk_file = varargin{2};
    pore_real_body_throat_Index = reshape(labelled_output, size_input);
    step6_vtkwrite(vtk_file, 'structured_points', 'Sample', pore_real_body_throat_Index);
end
waitbar(0.99,h,"All computations done, saving file");

% save results into output matfile
save (output_fname, 'num_pore_body', 'num_throat', 'labelled_output', 'pore_body_sizes', 'pore_throat_sizes', 'constriction_length', 'lsi', 'tortuosity', 'throat_wall_normalized_roughness', '-v7.3')
close (h);
%------------------------------------------------------------------- END CODE-------------------------------------------------------------------------------
