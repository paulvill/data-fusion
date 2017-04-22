function IM1 = rotate_image(IM2, R)
%ROTATE_IMAGE Rotate image using rotation matrix.
% IM1 = ROTATE_IMAGE(IM2, R) rotates the image im2 using rotation matrix R
% IM2 is the image to rotate
% 
% R is the 2x2 rotation matrix
% 
% IM1 is the rotated image

%%

% number of pixels
npixels = size(IM2, 1);

% calculate angle from rotation matrix
theta = atan2d(R(2,1), R(2,2));

% reference image
RI = imref2d([npixels npixels],[-0.5 0.5],[-0.5 0.5]);

% calculate corresponding affine transformation
A = affine2d([cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1]);

% transform image
IM1 = imwarp(IM2, RI, A, 'outputview', RI);

