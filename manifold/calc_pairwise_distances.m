% CALC_PAIRWISE_DISTANCES Calculate pairwise distances between images
%
% Usage
%    D = calc_pairwise_distances(images);
%
% Input
%    images: An array of images, indexed by the last dimension.
%
% Output
%    D: A square matrix of containing the pairwise Euclidean distances between
%       the images.

function D = calc_pairwise_distances(images)
    dim = ndims(images)-1;

    % Need to cast to double for calculations to work out if images are in
    % integer values.
    images = double(images);

    images = reshape(images, [], size(images, dim+1));

    n2 = sum(abs(images).^2, 1);

    D = bsxfun(@plus, n2', n2)-2*(images'*images);
    D = sqrt(max(0, D));
end

