function [inv_u] = ssl_estimate( W, D,indl, indu)
%SSL_ESTIMATE Summary of this function goes here
%   Detailed explanation goes here
	
	n = length(indl)+length(indu);
    indtot = [indl,indu];
    l = length(indl);
	u = n-l;
    T = zeros(length(indu)+length(indl));
    for i = 1:n
        for j = 1:n
            if j == indtot(i)
                T(i,j) = 1;
            end
        end
    end
    
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
