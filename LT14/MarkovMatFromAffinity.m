function [K] = MarkovMatFromAffinity( A )
%   MarkovMatFromAffinity generate a Markov matrix from an affinity matrix.
%
%       K_{ij} = A_{ij} / \sum_k(A_{kj})
%
%   Here, the sum of each COLUMN is required to be 1.
%
%   Input:
%    * A : A (n x n) affinity Matrix.
%
%   Output:
%    * K : A Markov Matrix.
%
%  (c) 2014  Roy R. Lederman and Ronen Talmon.

%
%   Sum of elements in each column.
%
    sum_d=sum(A);

%
%   Normalize the affinity matrix to be a Markov chain.
%
    K = A * diag(1./sum_d);

end

