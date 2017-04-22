function [MapEmbd, UWighted, d_t, singvals]=DiffusionMapsFromDistance(d,t,numNeighbors)
%
%   DiffusionMapsFromDistance computes diffusion distance and an embedding
%   from a matrix of pairwise distances using the Diffusion Maps algorithm.
%
%   For the details of the Diffusion Maps algorithm, see Appendix 1 of the
%   technical report.
%   The affinity matrix is constructed using a local scaling.
%
%   Input:
%    * d       : Symmetric matrix of pairwise distances.
%    * t       : Diffusion Time.
%    * numNeighbors : The number of neighbors to use in computing the local
%               scaling of the distances.
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
%                  A typical low dimensional representation of the embedding
%                  of the i-th point is  MapEmbd(i,[2,3]),
%                  where 2 and 3 are the first two non-trivial coordinates.
%    * UWeighted : The weighted eigenfunctions (in the columns), i.e.
%                  the approximate eigenfunctions of the
%                  Laplace-Beltrami operator.
%                  Note that the first column is the trivial eigenfunctions.
%    * d_t       : The diffusion distance at time t.
%    * singvals  : the singular values (eigenvalues) of the diffusion
%                   kernel.
%
%   The vectors are sorted in the columns of UWeighted, U,V in the order
%   of descending singular values (eigenvalues).
%
%   Typically, the useful output variables are MapEmbd and d_t.
%   The output variables UWeighted and singvals are typically used only in
%   the following context.
%   To compute the diffusion distance at a different time t',
%   compute the embedding at a different time:
%      M = UWighted * diag( svals.^t' )
%   and then, the diffusion distance at time t', between point i and j is
%   Euclidean distance between the corresponding rows of M:
%      d_{t'}(i,j) = || M(i,:) - M(j,:) ||_2
%
%  (c) 2014  Roy R. Lederman and Ronen Talmon.


    %
    % Compute the affinity matrix from the pairwise distances.
    %
    % A_ij = exp( d_{ij}^2 / ( sqrt{\epsilon_i} sqrt{\epsilon_j}),
    %   where \epsilon_i a scaling constant around the i-th point.
    %   \epsilon_i is taken to be the average of the square of the distances from
    %   the i-th point to its numNeighbors nearest neighbors.
    %
    [A,eps]=AffinityFromDistance( d, numNeighbors);
    %
    % Compute the normalized kernel for Diffusion Maps.
    %
    [K] = DiffusionKerFromAffinity( A );
    %
    % Spectral analysis of the kernel
    %
    [MapEmbd, UWighted, d_t, singvals, U,V] =DiffusionMapsFromKer( K , t );

end

