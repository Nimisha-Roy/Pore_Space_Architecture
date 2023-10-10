# Three-dimensional characterisation of pore space architecture in granular materials
This software segments the pore phase of a binary input microstructure into an interconnected network of bulky pore bodies and narrow pore throats, and computes various individual and collective properties of the pore space. Developed by Dr. Nimisha Roy.

## Citing this codebase
The program is written in MATLAB. It is developed by me (Nimisha Roy) as a part of my PhD work under the guidance of Dr. David Frost. The details of the software can be found in the citation referenced below.
If you use this codebase in your research or otherwise, please cite/reference the following paper:

Roy, N., Frost, J. D., & Roozbahani, M. M. (2023). Quantifying three-dimensional bodies and throats of particulate system pore space. Powder Technology, 415, 118160.

# Program Execution
Input: Matfile of a 3d array comprising solid phase and void phase labelled as 0 and 1 respectively
Give the input binary 3D structure, the software outputs the following properties of the structure:

1. 3D structure segmented into particle phase, individual pore bodies and individual pore throats
2. Sizes and shapes of individual pore bodies
3. Sizes and shapes of individual pore throats
4. Sizes of individual constrictions
5. Sharpness along the lengths of individual pore throats
6. Wall roughness of individual pore throats
7. Geometrical tortuosity of the structure in a given direction
8. Coordination number of pore bodies and pore throats

Copy the input structure matfile to the current folder and enter all input parameters in Master_Code.m. Run MAster_Code.m. All output parameters will be stored as a matfile in the same folder, along with an output segmented structure in vtk format for visualization. An example input structure called "test_input.mat" is present in the folder. The results from implementation of the code on the example input is saved as "test_output.mat" in the same folder.

## Input Parameters:
- fname (string) – filename of the input binary structure (should be a matfile with particle phase voxels designated as 0s and pore phase voxels designated as 1s).
- res (int) – resolution of the input structure
- output_fname (string) – output mat filename to store all results
- mpd (int) – mean particle diameter of input structure
- tortuous_axis (int, optional) - direction along which tortuosity is to be computed, if tortuosity is required as an output. 1 - along x axis, 2 - along y axis, 3 - along z axis.
- vtk_file (string, optional) - filename to save vtk output of labelled structure if needed.

## Output Parameters (stored in output matfile):
- num_pore_body (int) - number of pore bodies in the input structure
- num_throat (int) - number of pore throats in the input structure
- labelled_output (3D array) – labelled input structure, 0- particle phase, 1 to num_pore_body - Individual pore bodies, num_pore_body + 1 to num_pore_body + num_throat - Individual pore throats
- pore body sizes (2D array) - lengths of axes a,b,c of the ellipitic fit to the individual pore bodies
- pore throat sizes (2D array) - lengths of axes a,b,c of the ellipitic fit to the individual pore throats
- constriction sizes (1D array) - lengths of individual constriction surfaces
- lsi (1D array) - longitudinal sharpness indices of individual throats
- throat_wall_normalized_roughness (1D array) - normalized roughness of individual throat walls in percentage
- tortuosity (1D array) - tortuosity of paths along specified direction

## Output file for visualization:
- vtk_file.vtk - segmented and labelled output structure saved for visualization in Paraview

Please direct any questions or concerns to nroy9@gatech.edu.



