function image = image_normalize(image, rot_angle, symmetrize)

if nargin < 3 || isempty(symmetrize)
    symmetrize = false;
end

ndims0 = ndims(image);
if ndims0 == 3
    image = permute(image, [1 2 4 3]);
end

npixels = 100;

image = imresize(image, [npixels npixels]);
image = rotate_images(image, rot_angle);
for i = 1:size(image, 4)
    image(:,:,:,i) = adapthisteq(image(:,:,:,i), 'cliplimit',0.005);
    image(:,:,:,i) = crop_image(image(:,:,:,i), 0);
    if symmetrize
        image(:,:,:,i) = imlincomb(0.5, image(:,:,:,i), ...
            0.5, fliplr(image(:,:,:,i)));
        image(:,:,:,i) = adapthisteq(image(:,:,:,i));
    end
end

if ndims0 == 3
    image = permute(image, [1 2 4 3]);
end
