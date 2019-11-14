function [bool,init,time_init] = find_feasible_point(H,y,chi,xmin,xmax)
%======================================================================
% Solves P_I:            minimize    s  for (s in R, x in R^{N})
%                      subject to    -chi-s <= H*x-y <= chi+s
%                                    xmin <= x <= xmax
%                                    s >= 0
%
% Input----------------------------------------------------------------
%    H (matrix size qxN): observation operator (Radon transform)
%    y (vector length q): observed data
%    chi        (double): measurement uncertainty
%    xmin       (double): minimal pixel value
%    xmax       (double): maximal pixel value
%
% Output---------------------------------------------------------------
%    bool             (int): 1 if a feasible point was found, 0 else
%    init (vector length N): feasible geometry initialization  
%    time_init     (double): duration for solving the initialization pb 
%======================================================================

%%% problem P_I is re-written as:   minimize  c'*x 
%                                 subject to  A*x <= b
[M,N] = size(H);
A     =  [ -ones(M,1)            H;...
           -ones(M,1)           -H;...
           sparse(N,1)     speye(N);...
           sparse(N,1)    -speye(N);...
           -1               sparse(1,N)];
b     =  [y+chi;-y+chi;xmax.*ones(N,1);-xmin.*ones(N,1);0];
c     =  [1; sparse(N,1)];

% xc is a strictly feasible point for problem P_I 
xg  =  (xmax-xmin).*rand(N,1)+xmin;
s1  =  max(max(H*xg-y)-chi+1,0);
s2  =  max(max(y-H*xg)-chi+1,0);
xc  =  [max(s1,s2);xg];

%%% Newton barrier method starting at xc
mu = 10; % initialization of the barrier parameter
disp('------------------------------------------------------------------------')
disp('Start Newton barrier method to find a strictly feasible initial point')
disp('------------------------------------------------------------------------')
tic
x          = lp(A,b,c,xc,mu,1e-4,H,y,chi);
time_init  = toc;
init       = x(2:end);

%%% check if the geometry is strictly feasible
if xmin<min(init) && max(init)<xmax && max(abs(H*init-y))<chi
    bool=1;
    fprintf('\nFound a strictly feasible initial point after %.0f s.\n',time_init)
else
    bool=0;
    fprintf('\nDid not find a strictly feasible point.\n')
end
end