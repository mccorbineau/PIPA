function NormL = CalculNorme(L,x)
%========================================================================
% Computes the spectral norm of an operator using Algo. 12 on p.161 from 
%
% Pustelnik, Nelly. Méthodes proximales pour la résolution de problèmes 
% inverses: application à la tomographie par émission de positrons. 
% Diss. Université Paris-Est, 2010.
%
% Input -----------------------------------------------------------------
%    L     (function): operator
%    x (function arg): initial point 
%
% Output ----------------------------------------------------------------
%    NormL (double): spectral norm of L, ||L||
%========================================================================

%%% initialization
pnew    = 1 ;     % norm
nmax    = 20;     % maximal number of iterations
n       = 1 ;     % iteration number
epsilon = 0.01 ;  % required precision
cond    = 1;      % precision

%%% iterations
while (cond >= epsilon && n < nmax)
    xnew = L(x) ;
    p    = pnew ;
    pnew = norm(xnew) / norm(x) ;
    cond = abs(pnew-p)/ pnew ;
    x    = xnew;
    n    = n+1 ;
end
NormL = p;
end