function IMAGES2 = order_all_images(IMAGES, EMBED_COORD)
%ORDER_ALL_IMAGES Order a set of images by the value of embed_coord.
% IMAGES2 = ORDER_ALL_IMAGES(IMAGES, EMBED_COORD) orders the images stored
% in images by the value of embed_coord
% IMAGES is the array of images
% IMAGES is an npixels x npixels x nimages array (for 2D grayscale images)
%    npixels x npixels x nchannels x nimages array (for 2D color images)
%    npixels x npixels x nstack x nimages array (for 3D grayscale images)
%    npixels x npixels x nchannels x nstack x nimages array (for 3D color images)
% EMBED_COORD is an nimages-long vector
% IMAGES2 contains the ordered images; images2 is the same dimensions as
% the input array images

%%

% calulcate sorting indices
[~, I] = sort(EMBED_COORD);

% sort images
if ndims(IMAGES) == 3
    IMAGES2 = IMAGES(:,:,I);
elseif ndims(IMAGES) == 4
    IMAGES2 = IMAGES(:,:,:,I);
elseif ndims(IMAGES) == 5
    IMAGES2 = IMAGES(:,:,:,:,I);
end