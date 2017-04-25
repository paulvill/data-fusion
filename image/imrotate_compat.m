% IMROTATE Compatability layer for imrotate
%
% Usage
%    B = imrotate(A, angle, method, bbox)
%
% Input
%    A: Image to rotate.
%    angle: Rotation angle in degrees.
%    method: Interpolation method. One of 'nearest', 'bilinear', or 'bicubic'
%       (default 'nearest').
%    bbox: Bounding box type. Either 'loose' or 'crop' (default 'loose').
%
% Output
%    B: The image A rotated by angle.

function B = imrotate_compat(A, angle, method, bbox)
    if nargin < 3 || isempty(method)
        method = 'nearest';
    end

    if nargin < 4 || isempty(bbox)
        bbox = 'loose';
    end

    if isoctave
        B = imrotate(A, angle, method, bbox, 0);
    else
        B = imrotate(A, angle, method, bbox);
    end
end

