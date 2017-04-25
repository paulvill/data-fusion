% ZERO_PAD Zero-pad a signal
%
% Usage
%    x_zp = zero_pad(x, sz);
%
% Input
%    x: The signal to be zero-padded.
%    sz: The desired size of the output signal. Note that only the first
%       length(sz) dimensions will be zero-padded. The rest will be left
%       intact.
%
% Output
%    x_zp: The signal x, zero-padded along the first length(sz) dimensions
%       to have size sz. To recover the original signal, use the
%       extract_center function.

function x_zp = zero_pad(x, sz)
	[x, sz_roll] = unroll_dim(x, length(sz)+1);

	sz0 = size(x);
	sz0 = [sz0 ones(1, length(sz)-length(sz0))];
	sz0 = sz0(1:length(sz));

	if any(sz0 > sz)
		error('padding size must be larger than original signal');
	end

	x_zp = zeros([sz size(x, length(sz)+1)]);

	sub = cell(1, length(sz));
	for d = 1:length(sz)
		mid(d) = ceil((sz(d)+1)/2);
		ext1(d) = ceil((sz0(d)-1)/2);
		ext2(d) = floor((sz0(d)-1)/2);

		sub{d} = mid(d)-ext1(d):mid(d)+ext2(d);
	end
	sub{length(sz)+1} = 1:size(x, length(sz)+1);

	sub_grid = cell(size(sub));
	[sub_grid{:}] = ndgrid(sub{:});
	ind = sub2ind(size(x_zp), sub_grid{:});

	x_zp(ind) = x;

	x_zp = roll_dim(x_zp, sz_roll);
end

