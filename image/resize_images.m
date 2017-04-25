% RESIZE_IMAGES Resize a collection of images
%
% Usage
%    images = resize_images(images, npixels);
%
% Input
%    images: A set of images to be resized.
%    npixels: The desired size of the images.
%
% Output
%    images: The resized images.

function images_out = resize_images(images_in, npixels)
    if numel(npixels) == 1
        npixels = npixels*ones(1, 2);
    end

    sz_orig = size(images_in);

    sz_orig = [sz_orig ones(1, 4-length(sz_orig))];

    images_in = reshape(images_in, [sz_orig(1:2) sz_orig(3)*sz_orig(4)]);
    images_out = zeros([npixels sz_orig(3)*sz_orig(4)], class(images_in));

    for i = 1:size(images_in, 3)
        images_out(:,:,i) = imresize(images_in(:,:,i), npixels);
    end

    images_out = reshape(images_out, [npixels sz_orig(3) sz_orig(4)]);
end

