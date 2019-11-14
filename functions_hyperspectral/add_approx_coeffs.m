function [ Wx ] = add_approx_coeffs( P,J,nRow,nCol,approx,detail )
%=========================================================================
% Re-builds the full wavelet decomposition from separate approximation and 
% detail coefficients.
%
% Input -------------------------------------------------------------------
%    P         (int): number of endmembers (or materials, or channels)
%    J         (int): decomposition level   
%    nRow      (int): number of rows in the image
%    nCol      (int): number of columns
%    approx  (array): approximation coefficients of wavelet decomposition,
%                     size n_Row/(2^J) * n_Col/(2^J) * P 
%    detail  (array): detail coefficients of wavelet decomposition, 
%                     size nRow * nCol * P)
%
% Output ------------------------------------------------------------------
%    Wx (array nRow*nCol*P): full 2D wavelet decomposition
%==========================================================================

nRow_approx                            =  nRow/(2^J);
nCol_approx                            =  nCol/(2^J);
detail(1:nRow_approx,1:nCol_approx,:)  =  approx;
Wx                                     =  reshape(detail,nRow*nCol,P)';
end