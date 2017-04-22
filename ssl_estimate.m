function [inv_u] = ssl_estimate4( W, D, T,l)
%SSL_ESTIMATE Summary of this function goes here
%   Detailed explanation goes here
	
	n = size(W, 1);
	u = n-l;
    
    D_j = T*D*(T');
    W_j = T*W*(T');

	Dl = D_j(1:l,1:l);
	Du = D_j(l+1:l+u,l+1:l+u);

	Wll = W_j(1:l,1:l);
	Wlu = W_j(1:l,l+1:l+u);
	Wul = Wlu';
	Wuu = W_j(l+1:l+u,l+1:l+u);
    
    inv_u = (Du-Wuu)\Wul;

end
