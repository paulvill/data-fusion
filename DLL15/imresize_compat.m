%IMRESIZE_COMPAT Compatibility layer for imresize
% IM_OUT = IMRESIZE_COMPAT(IM_IN, ...) applies the imresize function to each
% channel of IM_IN, indexed along the third dimension. This is compatible with
% MATLAB's implementation of imresize, but Octave's version is restricted to
% only one or three channels. This compatibility wrapper allows for the use of
% any number of channels.
%
% For more information, see the documenation on imresize.

function im_out = imresize_compat(im_in, varargin)
    K = size(im_in, 3);

    im1 = imresize(im_in(:,:,1), varargin{:});

    im_out = zeros([size(im1, 1) size(im1, 2) K], class(im1));
    im_out(:,:,1) = im1;

    for k = 2:K
        im_out(:,:,k) = imresize(im_in(:,:,k), varargin{:});
    end
end

