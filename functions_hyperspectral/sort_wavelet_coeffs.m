function [approx,detail] = sort_wavelet_coeffs(Wx,nRow,nCol,P,J)
%==========================================================================
% Separates the detail coefficients from the approximation coefficients of 
% a given wavelet decomposition.
%
% Input -------------------------------------------------------------------
%    Wx    (matrix P*N): wavelet coefficients
%    nRow         (int): number of rows
%    nCol         (int): number of columns
%    P            (int): number of endmembers (or channels)
%    J            (int): decomposition level
%
% Output ------------------------------------------------------------------
%    approx (array): approximation wavelet coefficients
%    detail (array): detail wavelet coefficients
%==========================================================================

Wx           = reshape(Wx',nRow,nCol,P);
nRow_approx  = nRow/(2^J);
nCol_approx  = nCol/(2^J);

approx                                = Wx(1:nRow_approx,1:nCol_approx,:);
detail                                = Wx;
detail(1:nRow_approx,1:nCol_approx,:) = zeros(nRow_approx,nCol_approx,P);
end