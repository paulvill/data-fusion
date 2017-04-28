function [MapEmbd, UWighted, d_t, svals, U,V] = DiffusionMapsFromKer( K , t )
%   DiffusionMapsFromKer computes diffusion maps from a diffusion kernel.
%
%   Gives the embedding and diffusion distances for a given time t,
%   but also the eigenfunctions that can be used to construct the embedding
%   and distance for any other diffusion time.
%
%   For more details about the algorithm, see the algorithm for
%   approximating the  Laplace-Beltrami operator in
%      Lafon, S. , ``Diffusion Maps and Geometric Harmonics.''
%
%
%   Input:
%    * K  : An (n x n) diffusion kernel ( typically constructed
%           by DiffusionKerFromAffinity ).
%           Assumed to be symmetric, positive semi-definite.
%    * t  : The diffusion time
%
%   Output:
%    * MapEmbd  :  The diffusion maps embedding for time t.
%                  The i-th row MapEmbd(i,:) is the vector corresponding to
%                  the i-th point,
%                  The j-th column is the j-th coordinate of the
%                  embedding.
%                  Note, that the first column is a ``trivial''
%                  coordinate, and should be all ones. This coordinate can
%                  be ignored.
%    * UWeighted : The weighted eigenfunctions (in the columns), i.e.
%                  the approximate eigenfunctions of the
%                  Laplace-Beltrami operator.
%                  Note that the first column is the trivial eigenfunctions.
%    * d_t       : The diffusion distance at time t.
%    * seals     : the singular values (eigenvalues) of the diffusion
%                   kernel K.
%    * U,V       : The right and left singular vectors (eigenvectors) of K,
%                  such that
%
%   The vectors are sorted in the columns of UWeighted, U,V in the order
%   of descending singular values (eigenvalues).
%
%   To compute the diffusion distance at a different time t',
%   compute the embedding at a different time:
%      M = UWighted * diag( svals.^t' )
%   and then, the diffusion distance at time t', between point i and j is
%   Euclidean distance between the corresponding rows of M:
%      d_{t'}(i,j) = || M(i,:) - M(j,:) ||_2
%
%  (c) 2014  Roy R. Lederman and Ronen Talmon.

%
%   Eigenvalue decomposition of the kernel (assumed to be symmetric, positive
%   semidefinite)
%
    [U,S,V] = svd(K);
    svals=diag(S);

%
%   Compute the approximate eigenfunctions of the Laplace-Beltrami operator.
%
    UWighted = diag(1./U(:,1)) * U;

%
%   Compute the diffusion maps embedding at time t.
%
    MapEmbd = UWighted * diag( svals.^t );

%
%   Compute the diffusion distance at time t.
%
    d_t = distances(MapEmbd');

end

