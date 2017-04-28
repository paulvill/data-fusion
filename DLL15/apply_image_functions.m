function IMAGES = apply_image_functions(IMAGES_RAW, NPIXELS, DIM, CHANNEL_WEIGHT, CHANNEL_BLUR, CHANNEL_NORMALIZE, CHANNEL_MEAN_CENTER, RESIZE_IMAGE)
%APPLY_IMAGE_FUNCTIONS Apply relevant preprocessing to image set.
% IMAGES = APPLY_IMAGE_FUNCTIONS(IMAGES_RAW, DIM, CHANNEL_WEIGHT, CHANNEL_BLUR, CHANNEL_NORMALIZE, CHANNEL_MEAN_CENTER, RESIZE_IMAGE)
% applies several different preprocessing procedures to the images in IMAGES_RAW
% IMAGES_RAW is the array of images
% 
% IMAGES_RAW is an npixels x npixels x nimages array (for 2D grayscale images)
%    npixels x npixels x nchannels x nimages array (for 2D color images)
%    npixels x npixels x nstack x nimages array (for 3D grayscale images)
%    npixels x npixels x nchannels x nstack x nimages array (for 3D color images)
% 
% NPIXELS is resolution of the preprocessed images
% Images are assumed to be square, such that each preprocessed image will
% be NPIXELSxNPIXELS
%
% DIM is the dimension of each image
% DIM=2 means that the imgaes are each 2D
% DIM=3 means that the imgaes are 3-D z-stacks
% 
% CHANNEL_WEIGHT is a nchannel-long vector of weighting factors; each
% weighting factor should be between 0 and 1, where 0 means the channel
% does not contribute, and 1 means the channel is at its full strength
% 
% CHANNEL_BLUR is a nchannel-long vector of blurring factors; each
% weighting factor should be between 0 and 1, where 0 means no blurring and
% 1 means a blur radius which spans the entire image
% 
% CHANNEL_NORMALIZE is a nchannel-long Boolean vector which indicates which
% channels should be normalized; Contrast-limited Adaptive Histogram Equalization 
% via the adapthisteq function is applied to each channel where
% CHANNEL_NORMALIZE(i)=true
% 
% CHANNEL_MEAN_CENTER is a nchannel-long Boolean vector which indicates which
% channels should be used to mean-center the object in each image; 
% channels where CHANNEL_MEAN_CENTER(i)=true are used for edge-detection
% and mean-centering in the mean_center_image or mean_center_zstack
% functions
% 
% RESIZE_IMAGE is a Boolean variable which indicates whether the images
% should be resized; RESIZE_IMAGE=true means that each image is scale so
% that the object occupies 80% of the total image. 
% 
% IMAGES is the image set after processing


%% store necessary parameters

% allocate new images
IMAGES = imresize_compat(IMAGES_RAW, [NPIXELS NPIXELS]);

% number of images
nimages = size(IMAGES_RAW, ndims(IMAGES_RAW));

% store other parameters
% nstack =  number of images in one z-stack
% nchannels =  number of channels
if DIM == 3
    nstack = size(IMAGES_RAW, ndims(IMAGES_RAW)-1);
    if ndims(IMAGES_RAW) == 5
        nchannels = size(IMAGES_RAW, 3);
    else
        nchannels = 1;
    end
elseif ndims(IMAGES_RAW) == 4
    nchannels = size(IMAGES_RAW, 3);
else
    nchannels = 1;
end

% find the channels that should be used for mean-centering
channels = find(CHANNEL_MEAN_CENTER);

% create circular mask
% this will be applied to each image at the end (since the images are going
% to be rotated we will negelect the corners of the images)
[X, Y] = meshgrid(1:NPIXELS, 1:NPIXELS);
mask = ((X-(NPIXELS+1)/2).^2 + (Y-(NPIXELS+1)/2).^2 < ((NPIXELS+1)/2)^2);

%% apply image functions to each image

for i=1:nimages
    % 3D: extract relevant z-stack
    if DIM == 3
        if nchannels > 1
            im_tmp = IMAGES(:,:,:,:,i);
            for i2=1:nstack
                for j=1:nchannels
                    % apply image functions to each image
                    im_tmp(:,:,j, i2) = apply_image_functions_oneimage(im_tmp(:,:,j, i2), CHANNEL_NORMALIZE(j), CHANNEL_BLUR(j), CHANNEL_WEIGHT(j));
                end
            end
            % mean-center
            if ~isempty(channels)
                im_tmp = mean_center_zstack(im_tmp, channels, RESIZE_IMAGE);
            end
            % apply circular mask
            im_tmp = immultiply(im_tmp, repmat(mask, 1, 1, nchannels, nstack));
            
            IMAGES(:,:,:,:,i) = im_tmp;
        else
            im_tmp = IMAGES(:,:,:,i);
            for i2=1:nstack
                % apply image functions to each image
                im_tmp(:,:, i2) = apply_image_functions_oneimage(im_tmp(:,:,i2), CHANNEL_NORMALIZE, CHANNEL_BLUR, CHANNEL_WEIGHT);
            end
            % mean-center
            if CHANNEL_MEAN_CENTER == 1
                im_tmp = mean_center_zstack(im_tmp, 0, RESIZE_IMAGE);
            end
            % apply circular mask
            im_tmp = immultiply(im_tmp, repmat(mask, 1, 1, nstack));
            
            IMAGES(:,:,:,i) = im_tmp;
        end
    % 2D: extract relevant image
    else
        if nchannels > 1
            im_tmp = IMAGES(:,:,:,i);
            for j=1:nchannels
                % apply image functions to each image
                im_tmp(:,:,j) = apply_image_functions_oneimage(im_tmp(:,:,j), CHANNEL_NORMALIZE(j), CHANNEL_BLUR(j), CHANNEL_WEIGHT(j));
            end
            % mean-center
            if ~isempty(channels)
                im_tmp = mean_center_image(im_tmp, channels, RESIZE_IMAGE);
            end
            
            % apply circular mask
            im_tmp = immultiply(im_tmp, repmat(mask, 1, 1, nchannels));
            
            IMAGES(:,:,:,i) = im_tmp;
        else
            im_tmp = IMAGES(:,:,i);
            % apply image functions to each image
            im_tmp = apply_image_functions_oneimage(im_tmp, CHANNEL_NORMALIZE, CHANNEL_BLUR, CHANNEL_WEIGHT);
            % mean-center
            if CHANNEL_MEAN_CENTER == 1
                im_tmp = mean_center_image(im_tmp, 0, RESIZE_IMAGE);
            end
            
            % apply circular mask
            im_tmp = immultiply(im_tmp, mask);
            
            IMAGES(:,:,i) = im_tmp;
        end
    end
end

function IM2 = apply_image_functions_oneimage(IM1, NORMALIZE_SIGNAL, SIGNAL_BLUR_RADIUS, SIGNAL_SCALE)
%APPLY_IMAGE_FUNCTIONS_ONEIMAGE Applies relevant image functions to a
%single image
% IM2 = APPLY_IMAGE_FUNCTIONS_ONEIMAGE(IM1, NORMALIZE_SIGNAL, SIGNAL_BLUR_RADIUS, SIGNAL_SCALE)
% normalizes, blurs, and scales an image
% IM1 is the image to be altered
% NORMALIZE_SIGNAL is a Boolean; NORMALIZE_SIGNAL=true means Contrast-limited Adaptive Histogram Equalization (CLAHE)
% via the adapthisteq function is performed
%
% SIGNAL_BLUR_RADIUS is a blurring factors between 0 and 1, where 0 means no blurring and
% 1 means a blur radius which spans the entire image
%
% SIGNAL_SCALE is a weighting factors between 0 and 1, where 0 means the
% image is empty and 1 means the image is at its full strength

% normalize signal if required
if NORMALIZE_SIGNAL
    IM1 = adapthisteq(IM1);
end

% blur signal
npixels = size(IM1, 1);
blur_size = round(SIGNAL_BLUR_RADIUS*npixels);
if blur_size > 0
    filt = fspecial('disk', blur_size);
    IM1 = imfilter(IM1, filt, 'replicate');
end

% scale signal
IM2 = immultiply(IM1, SIGNAL_SCALE);

