function varargout = step6_rotation3dToEulerAngles(mat, varargin)
%This function extracts euler angles from a rotation matrix. Used in step6_ellipsoidfit of the algorithm
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2010-08-11,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.
%------------------------------------------------------------------- START CODE-------------------------------------------------------------------------------
p = inputParser;
validStrings = {'ZYX','ZYZ'};
addOptional(p,'convention','ZYX',@(x) any(validatestring(x,validStrings)));
parse(p,varargin{:});
convention=p.Results.convention;

% conversion from radians to degrees
k = 180 / pi;

switch convention
    case 'ZYX'
        % extract |cos(theta)|
        cy = hypot(mat(1,1), mat(2,1));
        % avoid dividing by 0
        if cy > 16*eps
            % normal case: theta <> 0
            phi   = k * atan2( mat(2,1), mat(1,1));
            theta = k * atan2(-mat(3,1), cy);
            psi   = k * atan2( mat(3,2), mat(3,3));
        else
            % 
            phi   = 0;
            theta = k * atan2(-mat(3,1), cy);
            psi   = k * atan2(-mat(2,3), mat(2,2));
        end
    case 'ZYZ'
        cy = hypot(mat(3,2), mat(3,1));
        if cy > 16*eps
            phi   = k * -atan2(mat(2,3), -mat(1,3));
            theta = k * -atan2(cy, mat(3,3));
            psi   = k * -atan2(mat(3,2), mat(3,1));
        else
            phi   = 0;
            theta = k * atan2(cy, mat(3,3));
            psi   = k * atan2(mat(2,1), mat(2,2));
        end
end

% format output arguments
if nargout <= 1
    % one array
    varargout{1} = [phi theta psi];
else
    % three separate arrays
    varargout = {phi, theta, psi};
end
%------------------------------------------------------------------- END CODE-------------------------------------------------------------------------------
