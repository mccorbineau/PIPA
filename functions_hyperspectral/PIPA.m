function [X,obj_vec,snr_global_vec,snr_mats_vec,time_vec] = PIPA(method,time_max,filename)
%==========================================================================
% Solves P_I:   minimize   0.5*||S*X-Y||^2_2 + reg*||W*X||_{1,detail} 
%                    for   (X in R^nEnd*nPix)
%             subject to   for every pixel p, sum_i(X_{i,p})<=1
%                          X >= 0
% where W is an orthogonal wavelet transform and ||W*.||_{1,detail} is the 
% l1 norm of the detail coefficients, using a proximal interior point algo-
% rithm.
%
% Kindly report any suggestions or corrections to
% mariecaroline.corbineau@gmail.com
%
% Input -------------------------------------------------------------------
%    method      (str): method, 'PIPA' (without variable metric) or
%                       'PIPA-VM' (with variable metric)
%    time_max    (int): maximal duration in seconds
%    filename    (str): name of the file to save the simulation results
%
% Output ------------------------------------------------------------------
%    X    (matrix nEnd*nPix): final abundance map
%    obj_vec        (vector): objective function value at every iteration
%    snr_global_vec (vector): global signal-to-noise ratio at every 
%                             iteration
%    snr_mats_vec   (vector): SNR for each material at every iteration
%    time_vec       (vector): time at every iteration
%
%==========================================================================
 
%%%%%%%% check the chosen method: PIPA or PIPA-VM (the latter includes a 
       % variable metric)
if strcmp(method,'PIPA')
    mode_metric = 0;
elseif strcmp(method,'PIPA-VM')
    mode_metric = 1;
else
    fprintf('This code does not handle method %s',method);
    return
end

%%%%%%%% problem parameters
[Y,S,nRow,nCol,Xtrue,materials] = load_data();
reg  = 0.01;        % regularization parameter
nPix = nRow*nCol;   % number of pixels
nEnd = size(S,2);   % number of materials or endmembers
C    = [speye(nEnd);-ones(1,nEnd)];    % constraint matrix
c    = [zeros(nEnd,1);1]*ones(1,nPix); % constraint vector    
StY  = S'*Y;
StS  = S'*S;                           % Hessian of the data-fitting term
if mode_metric==1
    % matrix involved in the Hessian of the barrier
    Cbig   = kron(speye(nPix),C);    
    % matrix involved in the Hessian of the data-fitting term  
    Qbig   = kron(speye(nPix),S'*S);   
    % minimal eigen value of the Hessian of the data-fitting term
    minStS = min(eig(S'*S));           
end

%%%%%%%% algorithm parameters
maxiter_backtrack = 1000;  % maximum number of iterations for the back-
                           % tracking search
theta_backtrack   = 0.80;  % granularity of backtracking search
delta_backtrack   = 0.99;  % Wolfe parameter of backtracking search
switch mode_metric
    case 1
        eps       = 1e5;   % precision for stopping criterion 
        gamma_max = 0.4;   % maximal step size
        mu        = 0.01;  % initial barrier parameter
    case 0
        eps       = 1e3;   % precision for stopping criterion 
        gamma_max = 1e-3;  % maximal step size
        mu        = 1;     % initial barrier parameter
end
zeta      = 1;     % used to ensure that stopping criteria decrease faster 
                   % than the barrier parameter
rho       = 1.5;   % geometric decrease for the barrier parameter    
precision = 1;     % precision for computing the proximity operator
VV        = [];    % used for warm restart when computing the prox

%%%%%%%% Wavelet operators
J        = 2;                             % resolution level
qmf      = MakeONFilter('Daubechies',4);  % wavelet 
[~,Jmax] = quadlength(rand(nRow,nCol));
% direct wavelet decomposition operator
Wav      = @(x) reshape(FWT2_PO(reshape(x,nRow,nCol),Jmax-J,qmf),nPix,1)'; 
% inverse wavelet decomposition operator
Wav_inv  = @(z) reshape(IWT2_PO(reshape(z,nRow,nCol),Jmax-J,qmf),nPix,1)';  

                             
%%%%%%%% useful functions
% gradient of smooth term + barrier
grad_phi = @(Z,CZ_c,mu_) -StY+StS*Z-mu_.*C'*(1./CZ_c);       
% proximity operator
if mode_metric==1
    prox     = @(Z,B,norm_B,step,prec,V) prox_metric(Z,B,norm_B,step,...
        prec,V,nRow,nCol,J,Wav,Wav_inv); 
else
    prox     = @(Z,step) prox_metric(Z,[],[],step,[],[],nRow,nCol,J,Wav,...
        Wav_inv);
end

%%%%%%%% initialization
counter       = 0;                    % used to save the results every hour
iter          = 1;
time          = 0;
X             = ones(nEnd,nPix)/(nEnd+1); % strictly feasible initial point
CX_c          = C*X+c;                    % constraints 
grad_phiX_old = grad_phi(X,CX_c,mu);      % gradient of the differentiable 
                                          % term + barrier
disp('------------------------------------------------------------------------')
fprintf('Start %s\n',method)
disp('------------------------------------------------------------------------')

while time<time_max
    if time-counter*3600>3600
        % save results every hour
        counter = counter + 1;
        save(filename,...
            'Y','S','reg','nRow','nCol','Xtrue','materials','X',...
            'obj_vec','snr_global_vec','snr_mats_vec','time_vec','time_max','method')
    end

    %%%%%%%% store variables to plot figures
    [snr_global_vec(iter),snr_mats_vec(:,iter)] = my_snr(X,Xtrue);
    obj_vec(iter)  = obj(X); % objective function 
    time_vec(iter) = time;   % running time                                        
    if (mod(iter-1,5)==0 && mode_metric==1) || (mod(iter-1,20)==0 && mode_metric==0)
        fprintf('iteration %d: f+g = %.3e | time = %.0f s | snr = %.2f dB\n',...
            iter-1,obj_vec(iter),time,snr_global_vec(iter))
    end
    tic        
    %%%%%%%% check stopping criterion to decrease the barrier parameter    
    if iter>1 && norm(RX(:))<eps*mu/zeta 
        mu        = mu/rho; % decrease barrier parameter
        %%% if there is a variable metric
        if mode_metric == 1 && mu<1e-6; rho= 1 + 1e-1; end
        if mode_metric == 1 && mu<4e-9; rho= 1 + 1e-2; end
        if mode_metric == 1 && mu<1e-12; rho= 1 + 1e-3; end
        %%% if there is no preconditioning
        if mode_metric == 0 && mu<5e-3; gamma_max = 2e-4; end
        if mode_metric == 0 && mu<5e-4; gamma_max = 8e-5; end
        
        % precision for computing the prox depends on barrier parameter
        precision = max(min(1,mu*1e3),1e-5); 
        % ensure that stopping criteria decrease faster than barrier parameter
        zeta      = zeta*(1+1e-5);           
    end
    
    if mode_metric==1
    %%%%%%%% build preconditioner
        Ssmall        =  mu./CX_c.^2;
        Sbig          =  diag(sparse(Ssmall(:)));  
        A             =  Qbig+Cbig'*Sbig*Cbig ; % preconditioner        
        A_1           =  @(u) A\u ;             % inverse of preconditioner
        A_1_grad_phiX =  reshape(A_1(grad_phiX_old(:)),nEnd,nPix);
        % norme of inverse of preconditioner
        N_A_1         =  1/(mu/max(max(CX_c(1:nEnd,:)))^2+minStS); 
    end
    
    %%% store current iterate
    Xold    = X;
    CXold_c = CX_c;
    
    %%%%%%%%% start backtracking
    gamma = gamma_max;
    for iter_backtrack=1:maxiter_backtrack
        if(reg>0)
            %%% forward-backward iteration
            if mode_metric==1
                [X,VV] = prox(Xold-gamma.*A_1_grad_phiX,A_1,N_A_1,...
                    gamma*reg,precision,VV);
            else
                X = prox(Xold-gamma.*grad_phiX_old,gamma*reg);
            end
        else
            if mode_metric==1
                X =  Xold-gamma.*A_1_grad_phiX;
            else
                X =  Xold-gamma.*grad_phiX_old;
            end
        end
        
        %%% check if the iterate is strictly feasible
        CX_c = C*X+c;
        if CX_c>0              
            %%% check the backtracking stopping criterion
            if mode_metric == 1
                upper_bound = delta_backtrack*((X(:)-Xold(:))'*A*(X(:)-Xold(:)))/(gamma*mu);
            else
                upper_bound = delta_backtrack*sum((X(:)-Xold(:)).^2)/(gamma*mu);
            end
            SXold_X   = S*(Xold-X);
            CX_CXold  = CX_c./CXold_c;
            if 0.5*sum(SXold_X(:).^2)/mu+sum(CX_CXold(:))-(nEnd+1)*nPix-sum(log(CX_CXold(:)))<=upper_bound 
                break
            end
        end
        if iter_backtrack>maxiter_backtrack
            fprintf('Primal variable backstracking did not converge\n')
            return
        end
        gamma = gamma*theta_backtrack; % decrease the stepsize
    end

    %%%%%%%% compute the inner loop stopping criterion
    grad_phiX = grad_phi(X,CX_c,mu);
    if mode_metric==1
        RX = A*(Xold(:) - X(:))./gamma - grad_phiX_old(:) + grad_phiX(:);
    else
        RX = (Xold(:) - X(:))./gamma - grad_phiX_old(:) + grad_phiX(:);
    end
    
    %%% update variables
    grad_phiX_old = grad_phiX;
    time          = time + toc;
    iter          = iter+1;
end

%%% store variables to plot figures 
% signal-to-noise ratio
[snr_global_vec(iter),snr_mats_vec(:,iter)] = my_snr(X,Xtrue); 
obj_vec(iter)  =  obj(X);  % objective function 
time_vec(iter) =  time;    % running time 
fprintf('iteration %d: f+g = %.3e | time = %.0f s | snr = %.2f dB\n',...
            iter-1,obj_vec(end),time,snr_global_vec(end))

%%% save results
save(filename,...
    'Y','S','reg','nRow','nCol','Xtrue','materials','X',...
    'obj_vec','snr_global_vec','snr_mats_vec','time_vec','time_max','method')
        
%%% recap  
disp('------------------------------------------------------------------------')
fprintf('Total duration for %d iterations of %s: %.1f s\n',iter,method,time)
fprintf('Final objective function value: %.3e\n',obj_vec(end))
fprintf('Final global SNR: %.2f dB\n',snr_global_vec(end))
fprintf('Final SNR for each material:\n')
for ii=1:length(materials)
fprintf('                             %s: %.2f dB\n',materials{ii},...
    snr_mats_vec(ii,end))
end
disp('------------------------------------------------------------------------')

function [OBJ] = obj(Z)
    %======================================================================
    % Computes the objective function value.
    % Input ---------------------------------------------------------------
    %    Z (matrix nEnd x nPix): abundance map
    % Output --------------------------------------------------------------
    %    OBJ (double): objective function value
    %======================================================================
    ZWt               = Wav_mult(Wav,Z);
    [~,detail_coeffs] = sort_wavelet_coeffs(ZWt,nRow,nCol,nEnd,J);
    SZ                = S*Z;
    OBJ = 0.5*norm(Y(:)-SZ(:))^2 + reg*(sum(abs(detail_coeffs(:))));
end
end