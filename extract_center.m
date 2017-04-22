% EXTRACT_CENTER Extract central part of signal
%
% Usage
%    x_cn = extract_center(x, sz);
%
% Input
%    x: The signal whose center should be extracted.
%    sz: The size of the central box to extract. Note that only the first
%       length(sz) dimensions of x will be indexed in this extraction. The
%       remaining dimensions will be untouched.
%
% Output
%    x_cn: The central box of size sz extracted from the first length(sz)
%       dimensions of x. To recover a signal of size x with the discarded
%       parts zeroed out, use the zero_pad function.

function x_cn = extract_center(x, sz)
	[x, sz_roll] = unroll_dim(x, length(sz)+1);

	% TODO: The similarities with zero_pad are significant. Needs to be
	% reduced.
	sz0 = size(x);
	sz0 = [sz0 ones(1, length(sz)-length(sz0))];
	sz0 = sz0(1:length(sz));

	if any(sz0 < sz)
		error('extraction box must be smaller than original signal')
	end

	x_cn = zeros([sz size(x, length(sz)+1)]);

	sub = cell(1, length(sz)+1);
	for d = 1:length(sz)
		mid(d) = ceil((sz0(d)+1)/2);
		ext1(d) = ceil((sz(d)-1)/2);
		ext2(d) = floor((sz(d)-1)/2);

		sub{d} = mid(d)-ext1(d):mid(d)+ext2(d);
	end
	sub{length(sz)+1} = 1:size(x, length(sz)+1);

	sub_grid = cell(size(sub));
	[sub_grid{:}] = ndgrid(sub{:});
	ind = sub2ind(size(x), sub_grid{:});

	x_cn = x(ind);

	x_cn = roll_dim(x_cn, sz_roll);
end

