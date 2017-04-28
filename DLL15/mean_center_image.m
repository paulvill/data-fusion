function IMAGE2 = mean_center_image(IMAGE, CHANNEL, RESIZE_IMAGE, edge_tol)
%MEAN_CENTER_IMAGE Translate an image so that it is (approximately)
%mean-centered in the frame.
% IMAGE2 = MEAN_CENTER_IMAGE(IMAGE, CHANNEL, resize_image, edge_tol)
% IMAGE is the image to be mean-centered 
% IMAGE is an npixels x npixels array (for grayscale images)
%    or an npixels x npixels x nchannels array (for color images)
% 
% CHANNEL is a vector of channels to use for the border detection for a color image
% (if the image is grayscale, the value of channel is ignored)
% 
% RESIZE_IMAGE is a Boolean variable
% RESIZE_IMAGE=true means that the object should be rescaled to occupy 80\%
% of the total image
% EDGE_TOL is the threshold used to detect edges in terms of how many edge
% pixels are needed in each row or column to qualify it as part of the
% object (default 0)
% 
% image2 is the cropped image

%% define relevant parameters

if nargin < 4 || isempty(edge_tol)
    edge_tol = 0;
end

% type of edge filer
detection_type = 'canny';

% number of pizels
npixels = size(IMAGE, 1);

%% compute edges in image

% extract relevant channels; store in im_tmp
if ndims(IMAGE) == 3
    im_tmp = IMAGE(:,:,CHANNEL);
    if length(CHANNEL) > 1
        im_tmp = uint8(mean(double(im_tmp), 3));
    end
else
    im_tmp = IMAGE;
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
    if ndims(IMAGE) == 3
        IMAGE2 = padarray(IMAGE, [nbuffer nbuffer 0]);
        IMAGE2 = imresize_compat(IMAGE2(idx3:(idx4+2*nbuffer), idx1:(idx2+2*nbuffer),:), [npixels npixels]);
    else
        IMAGE2 = padarray(IMAGE, [nbuffer nbuffer]);
        IMAGE2 = imresize_compat(IMAGE2(idx3:(idx4+2*nbuffer), idx1:(idx2+2*nbuffer)), [npixels npixels]);
    end
% otherwise, mean-center image
else
    RI = imref2d([npixels npixels],[1 npixels],[1 npixels]);
    tform = affine2d([1 0 0; 0 1 0; (1+npixels)/2-(idx1+idx2)/2 (1+npixels)/2-(idx3+idx4)/2 1]);
    IMAGE2 = imwarp(IMAGE, RI, tform, 'outputview', RI);
end

