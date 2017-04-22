function [EMBED_COORD, D2] = dm(W, EPS_SCALE, NCOMPS)
%DM Compute the diffusion maps embedding from the distance matrix W.
% [EMBED_COORD, D2] = DM(W, EPS_SCALE, NCOMPS) computes the diffusion maps
% embedding
% W is the matrix of pairwise distances
% 
% EPS_SCALE is the scaling of the median of the pairwise distances used for
% the kernel; EPS_SCALE=1 means the median of the pairwise distances are
% used for the diffusion maps kernel
% 
% NCOMPS is the number of embedding coordinates to calculate
% 
% EMBED_COORD contains the embedding coordinates (where the i^th column
% contains the i^th emebdding coordinate)
%
% D2 are the eigenvalues correspoinding to the embedding coordinates

%%

% epsilon for diffusion maps
eps = median(W(:)) * EPS_SCALE;

% number of eigenvectors to compute
neigs = NCOMPS + 1;

% calculate kernel of distances
W2 = exp(-(W/eps).^2);
W2 = diag(1./sum(W2)) * W2;

% compute eigenvectors
[V, D] = eigs(W2, neigs);

% sort by eigenvalue
[~, I] = sort(abs(diag(D)), 'descend');
V = V(:,I);
D = D(I,I);

EMBED_COORD = V(:, 2:end);
D2 = diag(D(2:end, 2:end));


