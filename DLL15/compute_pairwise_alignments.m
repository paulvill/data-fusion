function [R, W] = compute_pairwise_alignments(IMAGES, ANG_DIS)
%COMPUTE_PAIRWISE_ALIGNMENTS Compute the pairwise alignments (rotation
%matrices) and distances for use in a vector diffusion maps calculation.
% [R, W] = COMPUTE_PAIRWISE_ALIGNMENTS(images, ang_dis)
% compute the pairwise alignments and distances for a set of
% IMAGES, to be used in a vector diffusion maps calculation
% images is a set of images, where IMAGES(:,:,i) (for grayscale images) or
% IMAGES(:,:,:,i) (for color images) contains the i^th image in the dataset
% ANG_DIS is the angular discretication (in degreess) to use to compute the
% rotational pairwise alignments; smaller values of ang_dis will increase the
% required computational time, but could increase the accuracy
%
% R is a matrix containing the pairwise rotation matrices
% W is a matrix containing the minimum pairwise distances

%% store necessary parameters

% number of images
m = size(IMAGES, ndims(IMAGES));

% number of pixels
npixels = size(IMAGES, 1);

% dimension of rotation matrices
rot_dim = 2;

% angles to iterate over
theta_vec = 0:ANG_DIS:360;
theta_vec = theta_vec(1:end-1);
nrot = length(theta_vec);

% store rotations and distances
thetas = zeros(m);
R = zeros(rot_dim*m);
W = inf(m);

% reference image
RI = imref2d([npixels npixels],[-0.5 0.5],[-0.5 0.5]);


%% calculate pairwise alignments

h_waitbar = multi_waitbar(0,'Computing pairwise alignments...');

for i=1:nrot
    multi_waitbar(i/nrot, h_waitbar);
    theta = theta_vec(i);
    
    % calculate rotation transformation
    A = affine2d([cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1]);
    
    % apply rotation to images
    images_transformed = imwarp(IMAGES, RI, A, 'outputview', RI);
    
    % compute distances between transformed images and initial images
    dist_tmp = pdist2(reshape(double(IMAGES), [], m)', reshape(double(images_transformed), [], m)');
    
    % find images where this affine transformation is better than
    % all previous transforms
    idx = find(dist_tmp < W);
    
    % store new distances
    W(idx) = dist_tmp(idx);
    
    % store new rotation angles
    thetas(idx) = theta;
end

% compute rotation matrices from rotation angles
for i=1:m
    for j=1:m
        R(rot_dim*(i-1)+1:rot_dim*i, rot_dim*(j-1)+1:rot_dim*j) = calc_rot_matrix_2d(thetas(i, j));
    end
end

% symmetrize matrices
W = 0.5*(W + W');
R = 0.5*(R + R');

multi_waitbar(Inf, h_waitbar);


function R = calc_rot_matrix_2d(THETA)
% compute two-dimensional rotation matrix corresponding to a specific
% rotation
% THETA is the rotation angle
% R is the rotation matrix

R = [cosd(THETA) -sind(THETA);
    sind(THETA) cosd(THETA)];






