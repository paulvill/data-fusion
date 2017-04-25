%SSL_ESTIMATE Computes the matrix required to obtain the predicted sample
%values on unlabeled data points from labeled data points
%
% Usage 
%    [inv_u] = ssl_estimate( W, indl, indu)
%
% Input
%     W: affinity matrix containing the pairwise difference w_i,j between
%     all the data points
%     indl: indices of the labeled points
%     indu: indices of the unlabeled points
% 
% Output
%     inv_u: it is the matrix necessary to compute the predicted labels
%     from the unlabeled ones, f_u = inv_u*Y where Y correspond to the
%     labels

function [inv_u] = ssl_estimate( W, indl, indu)
	n = length(indl)+length(indu);
    indtot = [indl,indu];
    l = length(indl);
	u = n-l;
    
    % matrix transforming arbitrary distribution of labeled and unlabeled
    % data points into the classical ordered distribution of (un)labeled
    % data points, it simplifies the subsequent formulation of the solution
    T = zeros(length(indu)+length(indl));
    for i = 1:n
        for j = 1:n
            if j == indtot(i)
                T(i,j) = 1;
            end
        end
    end
    
    D = diag(sum(W, 2));
    
    D_j = T*D*(T');
    W_j = T*W*(T');

	Dl = D_j(1:l,1:l);
	Du = D_j(l+1:l+u,l+1:l+u);

	Wll = W_j(1:l,1:l);
	Wlu = W_j(1:l,l+1:l+u);
	Wul = Wlu';
	Wuu = W_j(l+1:l+u,l+1:l+u);
    
    % solution of the harmonic extension problem
    inv_u = (Du-Wuu)\Wul;

end
