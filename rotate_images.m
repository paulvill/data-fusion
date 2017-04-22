% ROTATE_IMAGES Rotate a collection of images
%
% Usage
%    images = rotate_images(images, angles);
%
% Input
%    images: The original images to rotate. This could either be a 3-dimensional
%       array consisting of a sequence of monochrome images, or a 4-dimensional
%       array made up of RGB images. In the former, the sequence is indexed
%       along the third dimension while in the latter, the fourth dimension is
%       used.
%    angles: An array of angles of same length as the number of images, specify-
%       ing the angles through which the images are to be rotated. If a single
%       angle is specified, all images are rotated through the same angle.
%
% Output
%    images: The rotated images.

function images = rotate_images(images, angles)
    ndims0 = ndims(images);
    if ndims0 == 3
        images = permute(images, [1 2 4 3]);
    end

    if numel(angles) == 1
        angles = angles*ones(size(images, 4), 1);
    end

    n_images = size(images, 4);

    for i = 1:n_images
        images(:,:,:,i) = imrotate_compat(images(:,:,:,i), angles(i), ...
            'nearest', 'crop');
    end

    if ndims0 == 3
        images = permute(images, [1 2 4 3]);
    end
end

