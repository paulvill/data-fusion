% LOCAL_NORMALIZE_IMAGES Locally normalize images by average value
% 
% Usage
%    images = local_normalize_images(images, sigma, epsilon);
% 
% Input
%    images: The images to be normalized in an N-by-N-by-M array.
%    sigma: The sigma parameter of the Gaussian kernel used to average the
%       images.
%    epsilon: The regularization parameter used when dividing by the local
%       average value to avoid singularity.
% 
% Output
%    images: The locally normalized images. At each point, the local average
%       avg is calculated using the Gaussian kernel of width sigma. The value
%       of the image at that point is then multiplied by y/(y^2+epsilon^2).

function images = local_normalize_images(images, sigma, epsilon)
    blur_ker = fspecial('gaussian',round(3*sigma), sigma);

    blur_nrm = @(im)(feval(@(x,y)(x.*y./(y.^2+epsilon.^2)), ...
        im, conv2(im, blur_ker, 'same')));

    images = matfun(blur_nrm, images, 3);
end

