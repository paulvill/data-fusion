function IMAGE_SET2 = register_all_images(IMAGE_SET, R_OPT)
%REGISTER_ALL_IMAGES Register all images using rotation matrices.
% IMAGE_SET2 = REGISTER_ALL_IMAGES(IMAGE_SET, R_OPT) registers all the
% images in IMAGE_SET using the rotations in R_opt calculated from vector
% diffusion maps
% IMAGE_SET is the array of images
% IMAGE_SET is an npixels x npixels x nimages array (for 2D grayscale images)
%    npixels x npixels x nchannels x nimages array (for 2D color images)
%    npixels x npixels x nstack x nimages array (for 3D grayscale images)
%    npixels x npixels x nchannels x nstack x nimages array (for 3D color images)
% 
% R_OPT is a 2*nimagesx2 matrix, whichcontains the optimal rotation matrices returned by the vdm
% function
% 
% IMAGE_SET2 is the images from IMAGE_SET, registered using the rotations in R_OPT

%%

% dimension of rotations
dim = size(R_OPT, 2);

% number of dimensions in image set
ndims_image_set = ndims(IMAGE_SET);

% number of images
nimages = size(IMAGE_SET, ndims_image_set);

if size(R_OPT, 1) ~= dim*nimages
    disp('ERROR: size of R_opt and image_set do not match');
    return
end

% allocate space for rotated images
IMAGE_SET2 = zeros(size(IMAGE_SET), 'like', IMAGE_SET);

% rotate each image and store
if ndims_image_set == 5
    for i=1:nimages
        R_tmp = R_OPT(dim*(i-1)+1:dim*i, :)';
        IMAGE_SET2(:, :, :, :, i) = rotate_image(IMAGE_SET(:,:, :,:, i), R_tmp);
    end
elseif ndims_image_set == 4
    for i=1:nimages
        R_tmp = R_OPT(dim*(i-1)+1:dim*i, :)';
        IMAGE_SET2(:, :, :, i) = rotate_image(IMAGE_SET(:,:, :,i), R_tmp);
    end
elseif ndims_image_set == 3
    for i=1:nimages
        R_tmp = R_OPT(dim*(i-1)+1:dim*i, :)';
        IMAGE_SET2(:, :, i) = rotate_image(IMAGE_SET(:,:,i), R_tmp);
    end
else
    disp('ERROR: invalid number of dimensions in image_set')
    return
end

