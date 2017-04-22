function [SuppressedMat] = SuppressDiagonal(Mat)
%   SuppressDiagonal erases the diagonal of a square matrix.
%
%   Sets Mat_ii=0 for all i.
%
%   Input:
%    * Mat : A (n x n) symmetric square matrix.
%
%   Output:
%    * SuppressedMat : The matrix Mat, with the diagonal suppressed.
%
%  (c) 2014  Roy R. Lederman and Ronen Talmon.

%
%   Validity test
%
    n=size(Mat,1);
    if (size(Mat,2)~=n) 
        error('Input must be a square matrix.');
    end

%
%   Suppressing the diagonal
%
    SuppressedMat = Mat - diag(diag(Mat));

end


