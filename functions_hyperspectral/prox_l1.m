function p = prox_l1(X,gamma) 
%==========================================================================
% Computes the proximity operator of gamma*the l_1 norm (soft thresholding)
% Solves the following problem (P): argmin_Y gamma*||Y||_1+0.5*||X-Y||^2
%
% Input -------------------------------------------------------------------
%    X      (array): point at which we want to compute the prox
%    gamma (double): multiplicative coefficient
%
% Output ------------------------------------------------------------------
%    p (array): proximity operator of the l1 norm at X, solution to (P)
%==========================================================================

p = (abs(X)>gamma).*(X-gamma*sign(X));
end