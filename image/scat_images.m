% SCAT_IMAGES Scatter images
%
% Usage
%    S_images = scat_images(images, scat_M, scat_J);
%
% Input
%    images: Images whose scattering transform is to be calculated. These are in
%       the form of an array of size N-by-N-by-K.
%    scat_M: The order of the scattering transform (default 1).
%    scat_J: The scale of the scattering transform. Specifically, the averaging
%       window size of the scattering transform is of width 2^scat_J (default
%       6).
%
% Output
%    S_images: The scattering transforms of the images in the form of an array
%       of size M-by-K, where M is the dimension of the scattering transform.

function S_images = scat_images(images, scat_M, scat_J)
    if nargin < 2 || isempty(scat_M)
        scat_M = 1;
    end

    if nargin < 3 || isempty(scat_J)
        scat_J = 6;
    end

    if ~exist('scat')
        error('Please run addpath_scatnet in scatnet toolbox');
    end

    ave_angles = false;
    pad_factor = 2;

    K = size(images, 3);

    % Prepare padding.
    orig_sz = size(images(:,:,1));
    padded_sz = pad_factor*orig_sz;

    % Set up scattering transform and scattering filter parameters.
    scat_opt = struct();
    scat_opt.M = scat_M;

    filt_opt = struct();
    filt_opt.J = scat_J;

    % Generate the wavelet filters of the transform.
    Wop = wavelet_factory_2d(padded_sz, filt_opt, scat_opt);

    % Apply scattering to the first image after zero-padding it.
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

    % Allocate matrix for all scattering coefficients.
    S_images = zeros([size(S_images1, 2) size(S_images1, 3) size(S_images1, 1) K]);

    S_images(:,:,:,1) = permute(S_images1, [2 3 1]);

    % Calculate for the remaining images.
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

    S_images = reshape(S_images, [], K);
end

function C = tensor_mult(A, B)
    sz_orig = size(B);
    B = reshape(B, sz_orig(1), prod(sz_orig(2:end)));
    C = A*B;
    C = reshape(C, [size(A, 1) sz_orig(2:end)]);
end

