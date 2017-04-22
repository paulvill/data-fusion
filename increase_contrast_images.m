% INCREASE_CONTRAST_IMAGES Increases contrast images through logistic function
%
% Usage
%    images = increase_contrast_images(images, a, b);
%
% Input
%    images: The images whose contrast is to be increased.
%    a, b: The scale and offset of the logistic function used.
%
% Output
%    images: The images filtered through the function
%       1/(1+exp(-a*(x-b))) .

function images = increase_contrast_images(images, a, b)
    images = 1./(1+exp(-a*(images-b)));
end

