function AVG_IMAGES = compute_average_trajectory(IMAGES, NSTAGES, WINDOW_SCALE)
%COMPUTE_AVERAGE_TRAJECTORY Compute average developmental trajectory for a
%set of (registrered and ordered) images
% AVG_IMAGES = COMPUTE_AVERAGE_TRAJECTORY(IMAGES, NSTAGES, WINDOW_SCALE)
% computes the average developmental trajectory
% IMAGES is a set of images, where images(:,:,:,i) is the i^th image
% NSTAGES is the number of images in the average trajectory
% WINDOW_SCALE the kernel scale for computing the average images
% AVG_IMAGES is the average trajectory, where AVG_IMAGES(:,:,:,i) is the
% i^th image in the average trajectory

%%

nimages = size(IMAGES, ndims(IMAGES));

avg_images_size = size(IMAGES);
avg_images_size(end) = NSTAGES;

% % check that images contains RGB images
% if ndims(IMAGES) ~= 4
%     disp('ERROR: average trajectory can only be computed over RGB images');
%     return
% end

AVG_IMAGES = zeros(avg_images_size, 'uint8');

% points at which to compute the average images, relative to the true
% images
frame_points = linspace(1, nimages, NSTAGES);

for i=1:NSTAGES
    
    % compute weights for averaging images, using a Gaussian kernel
    window_weights = exp(-((1:nimages)-frame_points(i)).^2/(WINDOW_SCALE/2)^2);
    window_weights = window_weights / sum(window_weights);
    
    % add images with relevant weights
    
    if ndims(IMAGES) == 5
        im_tmp = double(IMAGES(:,:,:,:,1))*window_weights(1);
        for j=2:nimages
            im_tmp = im_tmp + double(IMAGES(:,:,:,:,j))*window_weights(j);
        end
        AVG_IMAGES(:,:,:,:,i) = uint8(im_tmp);
    elseif ndims(IMAGES) ==4
        im_tmp = double(IMAGES(:,:,:,1))*window_weights(1);
        for j=2:nimages
            im_tmp = im_tmp + double(IMAGES(:,:,:,j))*window_weights(j);
        end
        AVG_IMAGES(:,:,:,i) = uint8(im_tmp);
    else
        im_tmp = double(IMAGES(:,:,1))*window_weights(1);
        for j=2:nimages
            im_tmp = im_tmp + double(IMAGES(:,:,j))*window_weights(j);
        end
        AVG_IMAGES(:,:,i) = uint8(im_tmp);
    end
%     im_tmp = immultiply(IMAGES(:,:,:,1), window_weights(1));
%     for j=2:nimages
%         im_tmp = imlincomb(1, im_tmp, window_weights(j), IMAGES(:,:,:,j));
%     end
    
    % store average image
%     AVG_IMAGES(:,:,:,i) = im_tmp;
    
end