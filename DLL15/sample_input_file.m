clear all
close all

%% Basic Algorithm
% The basic steps are as follows:
%
% # Read in images 
% # Preprocess the images (smooth and/or equalize channels, mean-center, etc.)
% # Run vector diffusion maps to register and order preprocessed images

%% Read Images
% The first step is to read in images.
% This is done using the |read_images| function 
% For two-dimensional images, it requires that all the images are stored in
% a single directory, all with a consistent prefix, and indexed with two
% digits starting from |01|. 
% For three-dimensional z-stack, it assumes each z-stack is stored in a
% single directory. 
% All of the z-stack directories are required to begin with a consistent prefix, and be indexed with two
% digits starting from |01|. 
% Similarly, all of the images are required to begin with a consistent prefix, and be indexed with two
% digits starting from |01| for each z-stack.

% directory where images are stored
image_dir = 'drosophila_fixed_images';

% prefix for each image
image_name = 'emb';

% image type/extension
image_ext = 'tif';

% no stack name or number of iamges in a stack are required for 2D images
stack_name = '';
nstack = 0;

% number of images
nimages = 120;

% dimension of images (dim=2 indicates standard 2D images, rather than
% z-stacks)
dim = 2;

% read in images
% images are stored in the variable images_raw
% nchannels stores the number of channels per image
[images_raw, nchannels] = read_images(image_dir, image_name, image_ext, ...
    stack_name, nimages, nstack, dim);

%% Show images
% An optional step is to now show the images using the |plot_images|
% function.

% plot the images
% images_raw are the images to plot 
% (returned from the read_images function)
% dim is the dimension of the images 
% (dim=2 indicates standard 2D images, rather than z-stacks)
plot_images(images_raw, dim)

%% Preprocess Images
% Now, we must preprocess the images using the |apply_image_functions| 
% function before registration and ordering, to
% remove any imaging and/or experimental artifacts.

% number of pixels
% images will be reduced to npixels x npixels resolution
npixels = 100;

% channel weights
% we scale the first (red) channel by half, and keep the second (green) and
% third (blue) channels at their input values
channel_weight = [0.5 1 1];

% channel blur
% we blur each of the channels by 5%
channel_blur = [0.05 0.05 0.05];

% channel normalization
% we normalize the first (red) channel using histogram equalization
% we do not normalize the second (green) or third (blue) channels
channel_normalize = [1 0 0];

% mean-center
% we use the first (red) channel to detect the edges of the object in order
% to mean center the object
channel_mean_center = [1 0 0];

% resize
% we choose to resize the images so all objects are (approximately) 
% the same size, to remove any variations due to size effects 
resize_image = true;

% we then apply these image functions of normalization, blurring,
% reweighting, and mean-centering
images = apply_image_functions(images_raw, npixels, dim, channel_weight, ...
    channel_blur, channel_normalize, channel_mean_center, resize_image);

% plot the images (optional)
plot_images(images, dim)

%% Calculate pairwise alignments
% We now need to calculate the angles needed to align _pairs_ of images

% angular discretization when computing pairwise aligments
% this means we search for pairwise aligmemnts over 10 degree increments
ang_dis = 10;

% compute the pairwise alignments 
% images are the preprocessed images
% ang_dis is the angular discretization
% R and W store the pairwise alignments and distances, respectively, for
% vector diffusion maps
[R, W] = compute_pairwise_alignments(images, ang_dis);

%% Apply vector diffusion maps
% We can now use vector diffusion maps to register and order the images. 

% ncomps is the number of components to compute
% we only compute 1 coordinate because we only need to order the images
% (i.e., sort by the first coordinate)
ncomps = 1;

% epsilon scale for diffusion maps kernel
% eps_scale = 0.25 means that the epsilon in the diffusion maps kernel is
% 1/4 of the median of the pairwise distances between data points
eps_scale = 0.25;

% vector diffusion maps calculates optimal rotations and embedding
% coordinate
[R_opt, embed_coord, D2] = vdm(R, W, eps_scale, ncomps);

% register images
images_registered = register_all_images(images, R_opt);

% order registered images by embedding coordinate
images_analyzed = order_all_images(images_registered, embed_coord);

% plot the images (optional)
plot_images(images_analyzed, dim)

%% Calculate average trajectory
% We can calculate an average trajectory from our set of (registered and
% ordered) images

% nsubimages is the desired number of images in the average trajectory
nsubimages = 10;

% avg_width is the (approximate) width of the averaging window used to
% compute each of the images in the average trajectory
avg_width = 4;

% compute the average trajectory
avg_images = compute_average_trajectory(images_analyzed, nsubimages, avg_width);

% plot the images (optional)
plot_images(avg_images, dim)
