function THETA = R_to_theta(R_OPT)
%R_TO_THETA Convert a set of rotation matrices to angles.
% theta = R_TO_THETA(R_OPT) converts the rotation matrices calculated from vector diffsion maps to rotation angles
% R_OPT is a 2nx2 matrix, where each 2x2 block corresponds to a rotation
% matrix
% THETA is a n-long vector containing the angles correponding to each
% rotation matrix in R_OPT

%% get relevant parameters

[m, dim] = size(R_OPT);

if dim ~= 2
    disp('ERROR: R_opt is not of the correct dimensions');
    return
end

% calculate number of rotations
n = m / dim;
if mod(n, 1) ~= 0
    disp('ERROR: R_opt is not of the correct dimensions');
    return
end

% convert rotation matrices to angles
THETA = zeros(n, 1);
for i=1:n
    R_tmp = R_OPT(dim*(i-1)+1:dim*i, :);
    THETA(i) = atan2d(R_tmp(2,1), R_tmp(1,1));
end


