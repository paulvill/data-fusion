function plot_images(IMAGES, IMAGE_DIM)
%PLOT_IMAGES Plot a set of images.
% PLOT_IMAGES(IMAGES, IMAGE_DIM) plots the images stored in images
% IMAGES is the array of images
% IMAGES is an npixels x npixels x nimages array (for 2D grayscale images)
%    npixels x npixels x nchannels x nimages array (for 2D color images)
%    npixels x npixels x nstack x nimages array (for 3D grayscale images)
%    npixels x npixels x nchannels x nstack x nimages array (for 3D color images)
% IMAGE_DIM is the dimension of each image
% IMAGE_DIM=2 means that the imgaes are each 2D
% IMAGE_DIM=3 means that the imgaes are 3-D z-stacks
% if IMAGE_DIM=3, then the maximum projection for each image is plotted

%%

% number of images
nimages = size(IMAGES, ndims(IMAGES));

% subplot dimensions
subplot_dim1 = round(sqrt(nimages));
subplot_dim2 = ceil(nimages/subplot_dim1);

% plot images
figure;
for i=1:nimages
%     subplot(subplot_dim1, subplot_dim2, i)
    make_subplot(subplot_dim1, subplot_dim2, i);
    % plot maximum projection if we have z-stacks
    if IMAGE_DIM == 3
        if ndims(IMAGES) == 4
            imshow(max_proj(IMAGES(:,:,:,i)))
        else
            imshow(max_proj(IMAGES(:,:,:,:,i)))
        end
    % otherwise, plot image
    else
        if ndims(IMAGES) == 3
            imshow(IMAGES(:,:,i))
        else
            imshow(IMAGES(:,:,:,i))
        end
    end
end


function IMAGE1 = max_proj(ZSTACK)
% calculcates maximum projection of ZSTACK

IMAGE1 = max(ZSTACK, [], ndims(ZSTACK));

function h = make_subplot(SUBPLOT_DIM1, SUBPLOT_DIM2, IDX)
% makes subplot with small margins (better for showing images)

% calculate space between subplots (10% of subplot)
subplot_space = min(0.1/SUBPLOT_DIM1, 0.1/SUBPLOT_DIM2);

% calculate bottom left corner of subplot
X_start = mod(IDX-1, SUBPLOT_DIM1) / SUBPLOT_DIM1 + subplot_space/2;
Y_start = 1 - ceil(IDX/SUBPLOT_DIM1)/SUBPLOT_DIM2 + subplot_space/2;

% return handle to subplot at specific location
h = subplot('position', [X_start Y_start 1/SUBPLOT_DIM1-subplot_space 1/SUBPLOT_DIM2-subplot_space]);