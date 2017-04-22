function [K2, K1, s1, s0] =DiffusionKerFromAffinity( AffinityMat )
%   DiffusionKerFromAffinity generate a diffusion kernel K2 for diffusion
%   maps, from a given symmetric affinity matrix A.
%
%   For more details about the algorithm, see the algorithm for
%   approximating the  Laplace-Beltrami operator in
%      Lafon, S. , ``Diffusion Maps and Geometric Harmonics.''
%
%       s0(j)   = \sum_k ( A_{kj} )
%       K1_{ij} = A_{ij} / ( s0(j) s0(i) )
%       s1(j)   = \sum_k ( K1_{kj} )
%       K2_{ij} = K1_{ij} / ( sqrt{s1(j)}  \sqrt{1(i)} )
%
%
%   Input:
%    * A  : A (n x n) affinity Matrix.
%
%   Output:
%    * K2 : The (n x n) diffusion kernel for diffusion maps.
%    * K1, s1, s0 are also given as part of the output, but they can be
%        ignored.
%
%  (c) 2014  Roy R. Lederman and Ronen Talmon.

%
%   Sum of elements in each column in the affinity matrix.
%
    s0=sum(AffinityMat);

%
%   Normalize the affinity matrix.
%
    K1 = diag(1./(s0) ) * AffinityMat * diag(1./(s0) );

%
%   Sum of elements in each column of K1.
%
    s1=sum(K1);

%
%   Normalize K1 to be a diffusion kernel.
%
    K2 = diag(1./sqrt(s1) ) * K1 * diag(1./sqrt(s1) );

end

