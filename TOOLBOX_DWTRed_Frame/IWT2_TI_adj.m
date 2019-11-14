function x = IWT2_TI_adj(tiwt,L,qmf)
% IWT2_TI -- Invert 2-d translation invariant wavelet transform
%  Usage
%    x = IWT2_TI(TIWT,L,qmf)
%  Inputs
%    TIWT     translation-invariant wavelet transform table, (3*(J-L)+1)*n by n
%    L        degree of coarsest scale
%    qmf      quadrature mirror filter
%  Outputs
%    x        2-d image reconstructed from translation-invariant transform TIWT
%
%  See Also
%    FWT2_TI, IWT2_TIMedian
%

[D1,n] = size(tiwt);
J = log2(n);
D = J-L;

%normalisation
for i=1:D-1
    tiwt((i-1)*3*n+1:i*3*n,:)=(1/2^(D-i))*tiwt((i-1)*3*n+1:i*3*n,:); 
end

x = IWT2_TI_adj_tight(tiwt,L,qmf);
%display('coucou')
%
% Copyright (c) 1995. David L. Donoho and Thomas P.Y. Yu
%


%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:40 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu
