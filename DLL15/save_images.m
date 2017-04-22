function save_images(IMAGES, DIM, IMAGE_DIR, IMAGE_NAME, IMAGE_EXT, STACK_NAME)
%SAVE_IMAGES Save images in specified directory and subdirectories.
% SAVE_IMAGES(IMAGES, DIM, IMAGE_DIR, IMAGE_NAME, IMAGE_EXT, STACK_NAME)
% saves the images stored in the array images 
% IMAGES is the array of images
% IMAGES is an npixels x npixels x nimages array (for 2D grayscale images)
%    npixels x npixels x nchannels x nimages array (for 2D color images)
%    npixels x npixels x nstack x nimages array (for 3D grayscale images)
%    npixels x npixels x nchannels x nstack x nimages array (for 3D color images)
% 
% DIM is the dimension of each image
% DIM=2 means that the imgaes are each 2D
% DIM=3 means that the imgaes are 3-D z-stacks
% 
% IMAGE_DIR is the directory where the images will be written
% 
% IMAGE_NAME gives the image prefix (for 2D images) or folder prefix (for
% 3D z-stacks) for the images
% IMAGE_EXT is the extension of the images (tif, jpg, etc.)
% 
% STACK_NAME is the image prefix for each of the images in a z-stack; if
% there are only 2D images, then stack_name is ignored

%%

% number of images
nimages = size(IMAGES, ndims(IMAGES));

% make image directory if it does not exist
if ~exist(IMAGE_DIR, 'dir')
    mkdir(IMAGE_DIR);
end

% loop through and write each image
h = multi_waitbar(0, 'Saving images...');
if DIM == 2
    if ndims(IMAGES) == 3
        for i=1:nimages
            multi_waitbar(i/nimages, h);
            
            % create file name of new image
            filename = sprintf('%s/%s%02d.%s', IMAGE_DIR, IMAGE_NAME, i, IMAGE_EXT);
            
            % write image, checking if file aleady exists
            check_val = write_single_image(IMAGES(:,:,i), filename);
            if check_val == -1
                multi_waitbar(Inf, h);
                return
            end
        end
    else
        for i=1:nimages
            multi_waitbar(i/nimages, h);
            
            % create file name of new image
            filename = sprintf('%s/%s%02d.%s', IMAGE_DIR, IMAGE_NAME, i, IMAGE_EXT);
            
            % write image, checking if file aleady exists
            check_val = write_single_image(IMAGES(:,:,:,i), filename);
            if check_val == -1
                multi_waitbar(Inf, h);
                return
            end
        end
    end
    
elseif DIM == 3
    nstack = size(IMAGES, ndims(IMAGES)-1);
    if ndims(IMAGES) == 4
        for i=1:nimages
            multi_waitbar(i/nimages, h);
            
            % make stack directory if it does not exist
            if ~exist(sprintf('%s/%s%02d', IMAGE_DIR, IMAGE_NAME, i), 'dir')
                mkdir(sprintf('%s/%s%02d', IMAGE_DIR, IMAGE_NAME, i));
            end
            for j=1:nstack
                
                % create file name of new image
                filename = sprintf('%s/%s%02d/%s%02d.%s', IMAGE_DIR, IMAGE_NAME, i, STACK_NAME, j, IMAGE_EXT);
                
                % write image, checking if file aleady exists
                check_val = write_single_image(IMAGES(:,:,j,i), filename);
                if check_val == -1
                    multi_waitbar(Inf, h);
                    return
                end
            end
        end
    else
        for i=1:nimages
            multi_waitbar(i/nimages, h);
            
            % make stack directory if it does not exist
            if ~exist(sprintf('%s/%s%02d', IMAGE_DIR, IMAGE_NAME, i), 'dir')
                mkdir(sprintf('%s/%s%02d', IMAGE_DIR, IMAGE_NAME, i));
            end
            for j=1:nstack
                
                % create file name of new image
                filename = sprintf('%s/%s%02d/%s%02d.%s', IMAGE_DIR, IMAGE_NAME, i, STACK_NAME, j, IMAGE_EXT);
                
                % write image, checking if file aleady exists
                check_val = write_single_image(IMAGES(:,:,:,j,i), filename);
                if check_val == -1
                    multi_waitbar(Inf, h);
                    return
                end
            end
        end
    end
end
multi_waitbar(h);

function CHECK_VAL = write_single_image(IM1, FILENAME)
% writes image IM1 to file FILENAME, checking if file exists already
% if file already exits, it asks the user if they want to overwrite the
% file
% CHECK_VAL=-1 means that the user would like to abort writing all image
% files

if exist(FILENAME, 'file')
    choice = questdlg(sprintf('%s already exists. Would you like to overwrite the image?', FILENAME), ...
        'File already exists', ...
        'Yes',...
        'No',...
        'Cancel writing files',...
        'No');
    % Handle response
    switch choice
        case 'Yes'
            CHECK_VAL = 0;
            imwrite(IM1,FILENAME);
        case 'No'
            CHECK_VAL = 0;
            return
        case 'Cancel writing files'
            CHECK_VAL = -1;
            return
    end
else
    CHECK_VAL = 0;
    imwrite(IM1,FILENAME);    
end
