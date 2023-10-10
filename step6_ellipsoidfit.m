function ell = step6_ellipsoidfit(points)
% This function finds the equivalent ellipsoid fit to a set of 3D points.
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
%------------------------------------------------------------------- START CODE-------------------------------------------------------------------------------
% number of points
n = size(points, 1);

% compute centroid
center = mean(points);

% compute the covariance matrix
covPts = cov(points)/n;

% perform a principal component analysis with 2 variables, to extract equivalent axes
[U, S] = svd(covPts);

% extract length of each semi axis
radii = sqrt(5) * sqrt(diag(S)*n)';

% sort axes from greater to lower
[radii, ind] = sort(radii, 'descend');

% format U to ensure first axis points to positive x direction
U = U(ind, :);
if U(1,1) < 0
    U = -U;
    % keep matrix determinant positive
    U(:,3) = -U(:,3);
end

% convert axes rotation matrix to Euler angles
angles = step6_rotation3dToEulerAngles(U);

% concatenate result to form an ellipsoid object
ell = [center, radii, angles];
%------------------------------------------------------------------- END CODE-------------------------------------------------------------------------------
