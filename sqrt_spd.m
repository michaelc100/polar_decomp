function B = sqrt_spd(A)
%SQRT_SPD   Matrix Square root of sym pos def
%   [B, its] = sqrt_spd(A) computes the 
%   square root of a symmetric positive definite A 
%   by doing a Cholesky decomposition and 
%   then performing a polar decomposition.
%   The matrix square root of A is given by the 
%   hermitian factor of this decomposition

R = chol(A);
[U, H, its] = poldec(R);
B = H;