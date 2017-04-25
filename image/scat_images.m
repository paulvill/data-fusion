% SCAT_IMAGES Scatter images
%
% Usage
%    S_images = scat_images(images, Wop, ave_angles, pad_factor);
%
% Input
%    images: Images whose scattering transform is to be calculated.
%    Wop: Wavelet operators given by wavelet_factory_2d.
%    ave_angles (optional): Determine whether to average along
%       rotation angles, creating rotation invariance (default false).
%    pad_factor: The factor with which the size of the images is to be
%       multiplied when zero-padding prior to computing the scattering
%       transform (default 2).
%
% Output
%    S_images: The scattering transforms of the images.

function S_images = scat_images(images, Wop, ave_angles, pad_factor)
    if nargin < 3 || isempty(ave_angles)
        ave_angles = false;
    end

    if nargin < 4 || isempty(pad_factor)
        pad_factor = 2;
    end

    if ~exist('scat')
        error('Please run addpath_scatnet in scatnet toolbox');
    end

    K = size(images, 3);

    % Prepare padding.
    orig_sz = [size(images, 1) size(images, 2)];
    padded_sz = pad_factor*orig_sz;

    [S_images1, meta] = ...
        format_scat(scat(zero_pad(images(:,:,1), padded_sz), Wop));

    % Unpad scattering coefficients.
    S_images1 = permute(S_images1, [2 3 1]);
    S_images1 = extract_center(S_images1, ceil(orig_sz/2^meta.resolution(1)));
    S_images1 = permute(S_images1, [3 1 2]);

    if ave_angles
        [~,~,I] = unique(meta.j', 'rows');
        sum_mat = sparse(I(:), [1:length(I)]', ones(length(I), 1));
        ave_mat = diag(sum(sum_mat, 2))^(-1)*sum_mat;
        S_images1 = tensor_mult(ave_mat, double(S_images1));
    end

    S_images = zeros([size(S_images1, 2) size(S_images1, 3) size(S_images1, 1) K]);

    S_images(:,:,:,1) = permute(S_images1, [2 3 1]);

    for k = 2:K
        S_imagesk = ...
            double(format_scat(scat(zero_pad(images(:,:,k), padded_sz), Wop)));
        S_imagesk = permute(S_imagesk, [2 3 1]);
        S_imagesk = extract_center(S_imagesk, ceil(orig_sz/2^meta.resolution(1)));
        S_imagesk = permute(S_imagesk, [3 1 2]);
        if ave_angles
            S_imagesk = tensor_mult(ave_mat, S_imagesk);
        end
        S_images(:,:,:,k) = permute(S_imagesk, [2 3 1]);
    end
end

function C = tensor_mult(A, B)
    sz_orig = size(B);
    B = reshape(B, sz_orig(1), prod(sz_orig(2:end)));
    C = A*B;
    C = reshape(C, [size(A, 1) sz_orig(2:end)]);
end

