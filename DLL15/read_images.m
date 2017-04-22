function [IMAGES, NCHANNELS] = read_images(IMAGE_DIR, IMAGE_NAME, IMAGE_EXT, STACK_NAME, NIMAGES, NSTACK, DIM)
%READ_IMAGES Read in images stored in specified directory and
%subdirectories.
% [IMAGES, NCHANNELS] = READ_IMAGES(IMAGE_DIR, IMAGE_NAME, IMAGE_EXT, STACK_NAME, NIMAGES, NSTACK, NPIXELS, DIM)
% reads in the images stored in image_dir
% IMAGE_DIR is the directory where the images are stored
% 
% IMAGE_NAME gives the image prefix (for 2D images) or folder prefix (for
% 3D z-stacks) for the images
% 
% IMAGE_EXT is the extension of the images (tif, jpg, etc.)
% 
% STACK_NAME is the image prefix for each of the images in a z-stack; if
% there are only 2D images, then stack_name is ignored
% 
% NIMAGES is the number of 2D images, or the number of 3D z-stacks
% 
% NSTACK is the number of images in a single z-stack; if
% there are only 2D images, then stack_name is ignored
% 
% DIM is the image dimension: dim=2 for standard 2D images, dim=3 for
% z-stacks
%
% IMAGES is the array of images
% IMAGES is an npixels x npixels x nimages array (for 2D grayscale images)
%    npixels x npixels x nchannels x nimages array (for 2D color images)
%    npixels x npixels x nstack x nimages array (for 3D grayscale images)
%    npixels x npixels x nchannels x nstack x nimages array (for 3D color images)
% NCHANNELS is the number of channels in the images

%%
try
    h = multi_waitbar(0, 'Reading images...');
    
    % read in first image to see if it is grayscale or color
    if DIM == 2
        i = 1;
        filename = sprintf('%s/%s%02d.%s', IMAGE_DIR, IMAGE_NAME, i, IMAGE_EXT);
    else
        i = 1;
        j = 1;
        filename = sprintf('%s/%s%02d/%s%02d.%s', IMAGE_DIR, IMAGE_NAME, i, STACK_NAME, j, IMAGE_EXT);
    end
    im_tmp = imread(filename);
    npixels = size(im_tmp, 1);
    if npixels ~= size(im_tmp, 2)
        % if the images are not square, we squarify them by adding some
        % zeros
        npixels_max = max(size(im_tmp, 1),size(im_tmp, 2));
        im_tmp = padarray(im_tmp, [npixels_max-size(im_tmp, 1) npixels_max-size(im_tmp, 2)]);
%         close(h);
%         msgbox('Images are not square');
%         return
    end
    if ndims(im_tmp) == 2
        NCHANNELS = 1;
    else
        NCHANNELS = size(im_tmp, 3);
    end
    
    % allocate space for images
    if DIM == 2
        if NCHANNELS == 1
            IMAGES = zeros(npixels, npixels, NIMAGES, 'uint8');
        else
            IMAGES = zeros(npixels, npixels, NCHANNELS, NIMAGES, 'uint8');
        end
    else
        if NCHANNELS == 1
            IMAGES = zeros(npixels, npixels, NSTACK, NIMAGES, 'uint8');
        else
            IMAGES = zeros(npixels, npixels, NCHANNELS, NSTACK, NIMAGES, 'uint8');
        end
    end
    
    % read in each image
    for i=1:NIMAGES
        multi_waitbar(i/NIMAGES, h);
        if DIM == 2
            % create filename of new image
            filename = sprintf('%s/%s%02d.%s', IMAGE_DIR, IMAGE_NAME, i, IMAGE_EXT);
            
            % read image
            im_tmp = imread(filename);
            
            % resize image
            im_tmp = imresize_compat(im_tmp, [npixels npixels]);
            
            % store image
            if NCHANNELS == 1
                IMAGES(:, :, i) = im_tmp;
            else
                IMAGES(:, :, :, i) = im_tmp;
            end
        else
            for j=1:NSTACK
                % create filename of new image
                filename = sprintf('%s/%s%02d/%s%02d.%s', IMAGE_DIR, IMAGE_NAME, i, STACK_NAME, j, IMAGE_EXT);
                
                % read image
                im_tmp = imread(filename);
                
                % resize image
                im_tmp = imresize_compat(im_tmp, [npixels npixels]);
                
                % store image
                if NCHANNELS == 1
                    IMAGES(:, :, j, i) = im_tmp;
                else
                    IMAGES(:, :, :, j, i) = im_tmp;
                end
            end
        end
    end
    multi_waitbar(Inf, h);
    
catch ME
    if strcmp(ME.identifier, 'MATLAB:imagesci:imread:fileDoesNotExist')
        close(h);
        msgbox(sprintf('%s does not exist', filename));
        rethrow(ME);
    end
end
