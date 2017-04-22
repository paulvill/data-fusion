function [A,eps]=AffinityFromDistanceGlobal( d, numNeighbors)
%   AffinityFromDistance Computes an affinity matrix
%   from given distances.
%
%   Computes    A_ij = exp( d_{ij}^2 / ( sqrt{\epsilon_i} sqrt{\epsilon_j}),
%   where \epsilon_i a scaling constant around the i-th point.
%   \epsilon_i is taken to be the average of the square of the distances from
%   the i-th point to its numNeighbors nearest neighbors.
%
%   Input:
%    * d : A (n x n) symmetric square matrix of distances.
%          d(i,j) is the distance between the i-th and j-th point.
%    * numNeighbors : The number of neighbors to use in computing the
%          scaling constant \epsilon.
%
%   Output:
%    * A : The (n x n) affinity matrix.
%    * eps : a vector of dimensionality n,
%            the local scaling \epsilon used in the computation.
%
%
%   See also: AffinityFromDistance_GlobalNorm
%
%  (c) 2014  Roy R. Lederman and Ronen Talmon.

    n=size(d,1);
    if (size(d,2)~=n)
        error('ERROR: d must be a square matrix.');
    end

    if (numNeighbors > (n-1))
        disp('Warning, setting numNeighbors to n-1.');
        numNeighbors = n-1;
    end

%
%   Compute the normalization coefficients.
%   For each point i, \epsilon_i is the average of the square of the distance
%   to the numNeighbors nearest neighbors of the point.
%
    [B,IX]=sort(d,2);
    eps=(mean(B(:,[2:numNeighbors+1])'.^2));

    eps = mean(eps)*ones(size(eps));

    if (~isempty(find(eps<=0)))
       error(['ERROR: eps=0, might be because the distance between some ' ...
          ' points is 0.']);
    end

%
%   Compute the affinity matrix
%
    A =  exp( - diag(1./sqrt(eps)) * (d.^2) * diag(1./sqrt(eps)) );

end

