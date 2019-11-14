function [ xt,xg,VV ] = prox_TV_metric(u,lambda,A_1,norm_A_1,VV,A,precision,L,Lt,N,TV)
%==========================================================================
% Computes the proximity operator of the total variation in a metric
% induced by a preconditioning matrix.
% Solves P_I:    minimize      0.5*||[xt;xg]-u||^2_A + lambda*TV(xg)  
%                for (xt and xg in R^N)
%
% Implementation of algorithm 2 in
% A. Chambolle and T. Pock. A first-order primal-dual algorithm for convex
% problems with applications to imaging. Journal of Mathematical Imaging 
% and Vision, Vol. 40, No. 1, pp 120-145, 2011.
%
% Input -------------------------------------------------------------------
%    u        (vector 2N): point at which the prox must be computed
%    lambda      (double): multiplicative coefficient 
%    A_1       (function): inverse of the preconditioner
%    norm_A_1    (double): norm of the inverse of the preconditioner
%    VV       (vector 2N): warm restart for the dual variable
%    A         (function): preconditioner 
%    precision   (double): precision for computing the prox
%    L         (function): computes the image horizontal and vertical 
%                          gradients 
%    Lt        (function): adjoint operator of L
%    N              (int): number of pixels
%    TV        (function): computes the total variation 
%
% Output ------------------------------------------------------------------
%    xt  (vector N): vectorized texture
%    xg  (vector N): vectorized geometry
%    VV (vector 2N): warm restart for dual variable
%==========================================================================

if isempty(VV)
    VV    = zeros(2*N,1);
    N_min = 100; % minimal number of iterations
else
    N_min=0; % minimal number of iterations
end
N_max = 1000; % maximal number of iterations

%%%%% algorithm parameters
tau   = (1e3/sqrt(8*norm_A_1));
gamma =  0.99;
sigma =  gamma/(tau*8*norm_A_1);

%%%%% initialization
ut   = u(1:N);
ug   = u(N+1:end);
xbar = u;
x    = u;
xt   = x(1:N);
xg   = x(N+1:end);

crit_ug = lambda*TV(ug);
% function to be minimized
crit    = lambda*TV(xg)+0.5*sum((x-u).*A(xt-ut,xg-ug)); 

%%%%% iterations
for i = 1:N_max
    crit_old = crit;
    
    temp     = VV+sigma*L(xbar);
    VV       = temp-sigma*prox_l2(temp./sigma,lambda/sigma,N);
    xold     = x;
    x        = (tau.*u+x-tau.*A_1(Lt(VV)))./(1+tau);
    theta    = 1/sqrt(1+2*gamma*tau);
    tau      = theta*tau;
    sigma    = sigma/theta;
    xbar     = x+theta.*(x-xold);  
    
    xt       = x(1:N);
    xg       = x(N+1:end);
    crit     = lambda*TV(xg)+0.5*sum((x-u).*A(xt-ut,xg-ug));

    % stopping criterion  
    residu=(crit-crit_old)/crit;
    if i>N_min && abs(residu) < precision && crit<crit_ug;break;end
    if i==N_max;disp('prox did not converge');end
end  
end


