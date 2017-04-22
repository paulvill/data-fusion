function image2 = crop_image(image, channel)
% function crops the black border from an image
% image: image to be cropped
% edge_tol: tolerance for the black border
% channel: channel to use for the border detection for an RGB image (if the image is grayscale, channel can be arbitrairy)
% image2: cropped image

edge_tol = 0;
detection_type = 'canny'; 

npixels = size(image, 1);

if ndims(image) == 3
    im_tmp = edge(image(:, :, channel), detection_type);
else
    im_tmp = edge(image, detection_type);
end

row_sum = sum(im_tmp);
idx1 = find(row_sum > edge_tol, 1, 'first');
idx2 = find(row_sum > edge_tol, 1, 'last');

col_sum = sum(im_tmp, 2);
idx3 = find(col_sum > edge_tol, 1, 'first');
idx4 = find(col_sum > edge_tol, 1, 'last');

if ndims(image) == 3
    image2 = imresize(image(idx3:idx4, idx1:idx2,:), [npixels npixels]);
elseif ndims(image) == 2
    image2 = imresize(image(idx3:idx4, idx1:idx2), [npixels npixels]);    
end