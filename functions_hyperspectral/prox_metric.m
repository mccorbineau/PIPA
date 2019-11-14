function [PX,v] = prox_metric(X,A_1,norm_A_1,gamma,precision,v,nRow,nCol,J,Wav,Wav_inv)
%==========================================================================
% Computes the proximity operator of gamma* the l1 norm of the detail 
% wavelet coefficients within a metric induced by a preconditioner.
%
% Implementation of the algorithm from
% P.-L. Combettes, D. Dung, and B.C. Vu. Proximity for sums of composite
% functions. Journal of Mathematical Analysis and Applications, vol. 380,
% no. 2, pp. 680-688, 2011.  
%
% Solves problem (P): argmin_Y 0.5*||X-Y||^2_A 
%                                    + gamma*||Wav*Y||_{1,detail}
%                   = argmin_Y 0.5*||A^{1/2}*X-A^{1/2}*Y||^2 
%                                    + gamma*||Wav*Y||_{1,detail}
%                   = A^{-1/2}* argmin_Y' 0.5*||A^{1/2}*X-Y'||^2 
%                                    + gamma*||Wav*A^{-1/2}*Y'||_{1,detail}
%
% Input -------------------------------------------------------------------
%    X     (matrix P*N): point at which the prox must be computed
%    A_1     (function): inverse of the preconditioner   
%    norm_A_1  (double): norm of the inverse of the preconditioning matrix
%    gamma     (double): coeff for the prox
%    precision (double): precision for solving (P)
%    v     (vector P*N): dual variable used for warm-restart
%    nRow         (int): number of rows
%    nCol         (int): number of columns
%    J            (int): decomposition level
%    Wav     (function): wavelet transform operator
%    Wav_inv (function): inverse wavelet transform operators
%
% Output ------------------------------------------------------------------
%    PX (matrix P*N): solution to (P)
%    v  (vector P*N): dual variable used for warm-restart
%==========================================================================

[P,N] = size(X); % number of materials and number of pixels

if(isempty(A_1)) 
    %%%% in case there is no preconditioning matrix
    WX               = Wav_mult(Wav,X);                                % wavelet transform
    [approx,detail]  = sort_wavelet_coeffs(WX,nRow,nCol,P,J);          % identify detail and approx wavelet coeffs
    detail           = prox_l1(detail,gamma);                          % apply the prox of the l1 norm on the detail coeffs
    prox_WX          = add_approx_coeffs(P,J,nRow,nCol,approx,detail); % re-built full wavelet coeffs
    PX               = Wav_inv_mult(Wav_inv,prox_WX);                  % go back to image space
else
    %%%%% start Dual Forward-Backward algorithm
    x = X(:); 
    if (isempty(v))
        
        %%%% if there is no warm-restart
        v = zeros(size(x));
    end
    
    %%%% parameters of DFB
    rho   =  1/(norm_A_1);
    delta =  min(1,rho) - 1e-8;
    step  =  2*rho-delta ;
    
    NbIt  = 1000; % maximal number of iterations
    pxold = x;    % store current iterate for the stopping criterion
    
    for i = 1:NbIt
        px = x - A_1(v) ;
        u  = v + step *px;
        v  = u - step .*reshape(prox_metric(reshape(u,P,N)/step,[],[],...
            gamma/step,[],[],nRow,nCol,J,Wav,Wav_inv),P*N,1); 
        if(norm(px - pxold) < precision && i>1); break; end
        pxold = px;
    end
    PX = reshape(px,P,N);
end

 