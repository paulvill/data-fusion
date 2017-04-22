function IMAGE_STACK2 = mean_center_zstack(IMAGE_STACK, CHANNEL, RESIZE_IMAGE, edge_tol)
%MEAN_CENTER_ZSTACK Coherently translate a set of z-stacks so that they are (approximately)
%mean-centered in the frame.
% image2 = MEAN_CENTER_ZSTACK(image_stack, channel, resize_image)
% IMAGE_STACK is the image stack to be mean-centered 
% IMAGE_STACK is npixels x npixels x nstack array (for grayscale images)
%    or npixels x npixels x nchannels x nstack array (for color images)
%
% CHANNEL is a vector of channels to use for the border detection for a color image
% (if the image is grayscale, the value of channel is ignored)
%
% RESIZE_IMAGE=true means that the object should be rescaled to occupy 80\%
% of the total image
% 
% IMAGE_STACK2 is the mean-centered z-stack

%% define relevant parameters

if nargin<4,
% threshold used to detect edges
edge_tol = 15;
end
% type of edge filer
detection_type = 'canny';

% number of pixels
npixels = size(IMAGE_STACK, 1);

%% compute edges in image

% extract relevant channels; store in im_tmp
im_tmp = uint8(mean(double(IMAGE_STACK), ndims(IMAGE_STACK)));
if ndims(im_tmp) == 3
    im_tmp = im_tmp(:,:,CHANNEL);
    if length(CHANNEL) > 1
        im_tmp = uint8(mean(double(im_tmp), 3));
    end
end

% adjust im_tmp
% im_tmp = imadjust(im_tmp);

% edge-detection
im_tmp = edge(im_tmp, detection_type);

% find indices where edges start and end
row_sum = sum(im_tmp);
idx1 = find(row_sum > edge_tol, 1, 'first');
idx2 = find(row_sum > edge_tol, 1, 'last');

col_sum = sum(im_tmp, 2);
idx3 = find(col_sum > edge_tol, 1, 'first');
idx4 = find(col_sum > edge_tol, 1, 'last');

% crop image to resize if RESIZE_IMAGE=true
if RESIZE_IMAGE
    nbuffer = round(0.1*npixels);
    if ndims(IMAGE_STACK) == 4
        IMAGE_STACKtemp = padarray(IMAGE_STACK, [nbuffer nbuffer 0 0]);
        IMAGE_STACK2 = IMAGE_STACK;
        for i = 1:size(IMAGE_STACK,4),
            for j = 1:3,
                IMAGE_STACK2(:,:,j,i) = imresize(IMAGE_STACKtemp(idx3:(idx4+2*nbuffer), idx1:(idx2+2*nbuffer),j,i),[npixels npixels]);
            end
        end
        %IMAGE_STACK2 = imresize_compat(IMAGE_STACK2(idx3:(idx4+2*nbuffer), idx1:(idx2+2*nbuffer),:,:), [npixels npixels]);
    else
        IMAGE_STACK2 = padarray(IMAGE_STACK, [nbuffer nbuffer 0]);
        IMAGE_STACK2 = imresize_compat(IMAGE_STACK2(idx3:(idx4+2*nbuffer), idx1:(idx2+2*nbuffer),:), [npixels npixels]);
    end
% otherwise, mean-center image
else
    RI = imref2d([npixels npixels],[1 npixels],[1 npixels]);
    tform = affine2d([1 0 0; 0 1 0; (1+npixels)/2-(idx1+idx2)/2 (1+npixels)/2-(idx3+idx4)/2 1]);
    IMAGE_STACK2 = imwarp(IMAGE_STACK, RI, tform, 'outputview', RI);
end
