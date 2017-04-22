function [R_OPT, EMBED_COORD, D2] = vdm(R, W, EPS_SCALE, NCOMPS)
%VDM Calculate the optimal rotations + embedding coordinates using vector
%diffusion maps.
% [R_OPT, EMBED_COORD, D2] = VDM(R, W, EPS_SCALE, NCOMPS) calculates the
% vector diffusion maps embedding coordinates and optimal rotations
% R is the matrix of pairwise rotation matrices (returned by
% compute_pairwise_alignments.m)
%
% W is the matrix of pairwise distances (returned by
% compute_pairwise_alignments.m)
%
% EPS_SCALE is the scaling of the median of the pairwise distances used for
% the kernel; eps_scale=1 means the median of the pairwise distances are
% used for the diffusion maps kernel
%
% NCOMPS is the number of embedding coordinates to calculate
% 
% R_OPT is a matrix containing the optimal rotations
%
% EMBED_COORD contains the embedding coordinates (where the i^th column
% contains the i^th emebdding coordinate)
% 
% D2 are the eigenvalue products correspoinding to the embedding
% coordinates

%% store relevant parameters

% dimension of rotations
dim = size(R,1) / size(W,1);

% number of data points
[n, ~] = size(W);

if mod(dim, 1) ~= 0
    disp('ERROR: sizes of R and W are not compatible')
    return
end

% calculate number of eigenvalues to compute
neigs = dim*(NCOMPS + 1);

% kernel scale for diffusion maps
eps = median(W(:)) * EPS_SCALE;

%% construct relevant matrices

% calculate kernel of distances
W2 = exp(-(W/eps).^2);
% make diaganol entries 0
for i=1:n
    W2(i, i) = 0;
end
% row-normalize
W2 = diag(1./sum(W2)) * W2;

% calculate entries of vdm matrix (rotation times kernel distances)
R2 = zeros(size(R));
for i=1:n
    for j=1:n
        R2(dim*(i-1)+1:dim*i,dim*(j-1)+1:dim*j) = W2(i,j) * R(dim*(i-1)+1:dim*i,dim*(j-1)+1:dim*j);
    end
end

%% calculate eigenvectors

% compute eigenvectors
[V, D] = eigs(R2, neigs);

% make sure matrix is real-valued (sometimes not due to round-off errors)
imag_tol = 1e-4;
count = 0;
count_max = 100;
while norm(imag(D)) < imag_tol && norm(imag(V)) > imag_tol && count < count_max
    [V, D] = eigs(R2, neigs);
    count = count + 1;
end    
if count == count_max || norm(imag(D)) > imag_tol
    disp('ERROR: imaginary eigenvectors')
end

% sort eigenvectors/eigenvalues
[~, ind] = sort(abs(diag(D)), 'descend');
V = V(:, ind);
D = D(ind, ind);
% normalize eigenvectors
for i=1:neigs
    V(:,i) = V(:,i) / norm(V(:,i));
end

%% calculate rotations

% if necessary, switch eigenvector sign to get proper rotation
if det(V(1:dim,1:dim)) < 0
    V(:,dim) = -V(:,dim);
end

% calculate optimal rotations using SVD
R_OPT = V(:,1:dim);
for i=1:n
    [u, s, v] = svd(R_OPT(dim*(i-1)+1:dim*i,:));
    R_OPT(dim*(i-1)+1:dim*i,:) = u * v';
end
R_OPT = R_OPT * R_OPT(1:dim,1:dim)';

%% calculate embedding coordinates from eigenvectors
EMBED_COORD = zeros(n, floor(neigs/dim)-1);
D2 = zeros(floor(neigs/dim)-1, 1);
for i=1:size(EMBED_COORD, 2);
    % find coordinate with maximum variance for each block of dim
    % eignevectors
    var_coord = 0;
    for j1=1:dim
        for j2=dim*i+1:dim*(i+1)
            embed_coord_tmp = sum(reshape(V(:,j1),dim, []).*reshape(V(:,j2),dim, []))';
            var_tmp = var(embed_coord_tmp);
            if var_tmp > var_coord
                var_coord = var_tmp;
                EMBED_COORD(:, i) = embed_coord_tmp;
                D2(i) = D(j1,j1) * D(j2,j2);
            end
        end
    end
end




