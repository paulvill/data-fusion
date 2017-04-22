function [d_t, Kt, K]=SimpleDiffusion(d,t,numNeighbors)
%
%   SimpleDiffusion computes a simplified diffusion distance,
%   and the corresponding diffusion kernel.
%   In this function, the normalization of the distances in the computation
%   of the affinity matrix is done using a vectors \epsilon
%   (local scale normalization) .
%   For alternatives, see  SimpleDiffusion_GlobalNorm
%
%   First, the diagonal-suppressed (``non-lazy'') affinity matrix B
%   is computed from the distance matrix d,
%        B_{ij} = exp( d_{ij}^2 / (\sqrt{\epsilon_i) \sqrt{\epsilon_j}) )
%             where i \ne j,
%       and B_{ii} = 0.
%   Then, the affinity matrix is normalized to obtain the Markov kernel K,
%       K_{ij} = A_{ij} / \sum_k(A_{kj}) ,
%     note that here, the sum of each COLUMN is required to be 1.
%
%   Next, we compute t steps of the Markov chain,
%     Kt = K^t.
%   Finally, the diffusion distance at time t is the Euclidean distance
%   between the corresponding columns of Kt.
%
%
%   Input:
%    * d  : An (n x n) distance matrix.
%           Assumed to be symmetric, with non-negative elements.
%    * t  : The diffusion time
%
%   Output:
%    * d_t : The diffusion distance at time t.
%    * Kt   : The Markov matrix for the diffusion, after t steps.
%    * K   : The Markov matrix for the diffusion.
%        note that here, the sum of each COLUMN is required to be 1.
%
%  (c) 2014  Roy R. Lederman and Ronen Talmon.

%
%   Compute the affinity matrix, and set the diagonal to 0.
%
    [A,eps1]=AffinityFromDistance( d, numNeighbors);
    [B] = SuppressDiagonal(A);
    %B=A;     % uncomment to override the suppression of the diagonal.

%
%   Compute the Markov Kernel by normalizing the affinity matrix.
%
    K = MarkovMatFromAffinity( B );

%
%   t steps of the diffusion.
%
    Kt = K^t;

%
%   Compute the pairwise diffusion distance.
%

    d_t = pdist2(Kt',Kt','euclidean');

end

