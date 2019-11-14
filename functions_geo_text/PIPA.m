function [x,obj_vec,snr_vec,time_vec] = PIPA(sample_name,time_max,filename)
%======================================================================================================
% Solves P_I:      minimize   0.5*||F(xt)||^2_2 + weightTV*TV(xg)  
%                       for   (xt and xg in R^N)
%                subject to   -chi <= H*(xt+xg)-y <= chi
%                             xmin <= xt+xg <= xmax
%                             -xmax/3 <= xt <= xmax/3
% where TV is the isotropic total variation with Dirichlet conditions, F is an edge detection operator 
% derived form the Laplacian, using a proximal interior point algorithm.
%
% Kindly report any suggestions or corrections to
% mariecaroline.corbineau@gmail.com
%
% Input -----------------------------------------------------------------------------------------------
%    sample_name (str): sample name, 'agaricus' or 'glass'
%    time_max    (int): maximal total duration (including initialization) in seconds  
%    filename    (str): name of the file to save the simulation results
%
% Output ----------------------------------------------------------------------------------------------
%    x          (vector length 2N): vectorized geometry-texture decomposition, x(1:N)=texture, 
%                                   x(N+1:end)=geometry 
%    obj_vec              (vector): objective function value for every iteration
%    snr_vec              (vector): signal-to-noise ratio for every iteration
%    X_Xinf_vec           (vector): distance to solution for every iteration
%    time_vec             (vector): time for every iteration (includes duration for initialization)
%========================================================================================================

%%%%% Load data and feasible initialization point   
[x_true,y,chi,H,weightTV] = load_data(sample_name);
%   x_true (matrix size mxn): ground-truth (used only to compute the signal-to-noise ratio)
%   y      (vector length q): observed data
%   chi             (double): measurement uncertainty
%   H      (matrix size qxN): observation operator (Radon transform)
%   weightTV        (double): regularization parameter 

%%%%%%%%%%%% algorithm hyperparameters 
maxiter_backtrack  = 100;   % maximum number of iterations for the inner loop
mu                 = 1e-3;  % initialization of barrier parameter
precision          = 5e-6;  % precision for computing the proximity operator of TV
gamma_max          = 1.2;   % maximal stepsize
theta_backtrack    = 0.8;   % granularity of backtracking search
delta_backtrack    = 0.99;  % Wolfe parameter of backtracking search     
eps                = 1e4;   % accuracy for stopping criterion 
rho                = 1.1;   % constant for geometric decrease of barrier parameter 
zeta               = 1.2;   % constant for geometric decrease of accuracy
VV                 = [];    % used for warm-restart in the accelerated Chambolle-Pock algorithm

%%%%%%%%%%%% problem parameters
xmin    = 0;                   % minimal pixel value
xmax    = 1;                   % maximal pixel value
[m,n]   = size(x_true);        % image size
N       = size(H,2);           % CT operator size
Htildet = (H.*(H*ones(N,1)))'; % used to compute the preconditioning matrix

%%%%%%%%%%%% creation of the edge detection operator
% the FFT is used to compute efficiently the convolution
Laplacian  = [0 1 0;1 -4 1;0 1 0];
W          = padarray(Laplacian,[m-3+1,n-3+1]/2,'pre');
W          = padarray(W,[m-3-1,n-3-1]/2,'post');
W          = ones(m,n)-fft2(fftshift(W));                                  % edge detection operator
F          = @(x) reshape(W.*fft2(reshape(x,m,n)),N,1)./m;                 % convolution operator m=n
F_adj      = @(z) reshape(real(ifft2(W.*reshape(z,m,n))),N,1).*m;          % adjoint
F_adj_F_1  = @(x) reshape(real(ifft2(fft2(reshape(x,m,n))./(W.^2))),N,1);  % inverse of F_adj*F

%%%%%%%%%%%% definition of useful functions
smooth     = @(xt)     0.5*sum(abs(F(xt)).^2);           % smooth regularization  
obj        = @(xt,xg)  smooth(xt)+weightTV*calcTV(xg);   % objective function
my_snr     = @(xt,xg) -20*log10(norm(x_true(:)-xt-xg)/norm(x_true(:))); % signal-to-noise ratio
TV         = @(x)      calcTV(x);                        % total variation
L          = @(x)      LinOp(x);                         % gradient operator (vertical, horizontal)
Lt         = @(x)      LinOpT(x);                        % transpose of the gradient operator 
% proximity operator of TV
prox_TV = @(u,lambda,A_1,norm_A_1,VV,A,precision) ...
    prox_TV_metric(u,lambda,A_1,norm_A_1,VV,A,precision,L,Lt,N,TV); 

%%%%%%%%%%%% Find an feasible initial point
xt                  = zeros(N,1);                        % texture
[bool,xg,time_init] = find_feasible_point(H,y,chi,xmin,xmax); % initialization of geometry
%   xg     (vector length N): geometry component initialization    
%   time_init       (double): duration for solving the initialization problem 
%   bool               (int): 1 if a feasible initial point is found, 0 else

if bool==1
    [C1,C2,C3,C4,C5,C6,Cfull] = constraints(xt,xg); % constraints at initial point
    temp_old                  = mu.*(1./C1-1./C2+1./C5-1./C6+H'*(1./C3-1./C4));
    F_adj_F_xt_old            = F_adj(F(xt));
    gamma                     = 0;                  % stepsize
    iter                      = 1;
    time                      = time_init;
    counter                   = 0;                  % used to save the results every hour

    %%%%%%%%%%%% start PIPA
    disp('------------------------------------------------------------------------')
    disp('Start PIPA-VM')
    disp('------------------------------------------------------------------------')

    while time < time_max
        if time-counter*3600>3600
            % save results every hour
            x           = [xt;xg];
            counter     = counter + 1;
            save(filename,...
                'sample_name','time_max','x_true','y','chi','H','weightTV','time_init',...
                'obj_vec','snr_vec','time_vec','x','x_true')
        end
        %%%%%% store variable to plot figures
        obj_vec(iter)  = obj(xt,xg);                                      % objective function
        snr_vec(iter)  = my_snr(xt,xg);                                   % signal-to-noise ratio
        time_vec(iter) = time;                                            % time
        
        if (mod(iter-1,5)==0 && iter-1<20)||(mod(iter-1,10)==0 && iter-1<40)||mod(iter-1,100)==0 
            fprintf('iteration %d: f+g = %.3e | time = %.0f s | snr = %.2f dB\n',iter-1, obj_vec(iter),time,snr_vec(iter)); 
        end

        tic;
        %%%%%%%%%%%% check stopping criterion for Algorithm 1 
        if iter>1 && norm(RX(:))<eps*mu/zeta 
            mu        = mu/rho;                 % decrease barrier parameter
            zeta      = zeta*1.00001;           % accuracy must decrease faster than mu
            if mu>2e-4
                precision = mu/800; % precision for computing the prox of TV
            elseif mu>1.3e-4
                precision = 1e-7;
            else
                precision = 1e-8;
            end
        end  

        %%%%%%%%%%% build preconditioning matrix
        D1   = mu.*(1./(C1.^2)+1./(C2.^2)+1./(C5.^2)+1./(C6.^2));
        D2   = mu.*(1./(C3.^2)+1./(C4.^2));
        G    = D1 + Htildet*D2;
        A    = @(ut,ug) [F_adj(F(ut))+G.*(ut+ug); G.*(ut+ug)]; % preconditioner
        A_1  = @(x)     [F_adj_F_1(x(1:N)-x(N+1:end)); ...     % inverse of preconditioning operator
                             F_adj_F_1(x(N+1:end)-x(1:N)) + x(N+1:end)./G]; 
        % compute the norm of the inverse of the preconditioner
        norm_A_1 = CalculNorme(A_1,[xt;xg]);  
        
        %%%%%%%%%%% update iterates
        xt_old     = xt;    % texture
        xg_old     = xg;    % geometry
        Cfull_old  = Cfull; % constraints

        %%%%%%%%%%% start backtracking
        gamma = gamma_max; % stepsize
        for iter1=1:maxiter_backtrack  
            %%%%% forward-backward iteration
            x          = [xt_old;xg_old]-gamma.*[xt_old;-xt_old-temp_old./G];
            [xt,xg,VV] = prox_TV(x,gamma*weightTV,A_1,norm_A_1,VV,A,precision);
            [C1,C2,C3,C4,C5,C6,Cfull] = constraints(xt,xg);
             
             %%%%% check if the iterate is feasible
             if Cfull>0  
                 A_x_old_x   = A(xt_old-xt,xg_old-xg);
                 upper_bound = sum([xt_old-xt;xg_old-xg].*A_x_old_x)*delta_backtrack/gamma;
                 criterion   = smooth(xt)-smooth(xt_old)-mu*sum(log(Cfull./Cfull_old))+...
                     sum((xt+xg-xt_old-xg_old).*temp_old)-sum((xt-xt_old).*F_adj_F_xt_old);
       
                 %%%%% check if (5) is satisfied
                 if criterion <= upper_bound ; break; end 
             end
             
             if iter1==maxiter_backtrack ; disp('Backtracking did not converge'); x=[xt;xg]; return; end
             gamma = gamma*theta_backtrack; % decrease gamma
        end

        %%%%%%%%%%% compute intermediate variables for the stopping criterion
        temp       = mu.*(1./C1-1./C2+1./C5-1./C6+H'*(1./C3-1./C4));
        F_adj_F_xt = F_adj(F(xt));
        RX         = A_x_old_x./gamma + [F_adj_F_xt-F_adj_F_xt_old-temp+temp_old;-temp+temp_old];
         
        %%%%%%%%%%% update variables, time, and iteration number
        F_adj_F_xt_old = F_adj_F_xt;
        temp_old       = temp;
        time           = time + toc;  
        iter           = iter + 1;
    end
    x = [xt;xg];
    %%% store results
    obj_vec(iter)  =  obj(xt,xg);                                    % objective function
    snr_vec(iter)  =  my_snr(xt,xg);                                 % signal-to-noise ratio
    time_vec(iter) =  time;                                          % time
    %%% recap
    disp('------------------------------------------------------------------------')
    fprintf('Duration for initialization: %.1f s\n', time_init)
    fprintf('Total duration after initialization and %d iterations of PIPA: %.0f s\n', iter, time)
    fprintf('Final objective function value: %.3e\n', obj_vec(end))
    fprintf('Final SNR: %.2f dB\n', snr_vec(end))
    disp('------------------------------------------------------------------------')
    
    %%%%%%%%%%% save results
    save(filename,...
            'sample_name','time_max','x_true','y','chi','H','weightTV','time_init',...
            'obj_vec','snr_vec','time_vec','x','x_true')   
end
function tv = calcTV(x)
   % Computes the isotropic total variation of an image using Dirichlet boundary conditions.
   % Inputs -------
   %    x (vector length N): vectorized 2D image of size mxn
   % Outputs ------
   %    tv (double): isotropic total variation
   X        = reshape(x,m,n);                   
   U        = X;
   U(:,2:n) = X(:,2:n)-X(:,1:(n-1));    
   V        = X;
   V(2:m,:) = X(2:m,:)-X(1:(m-1),:);    
   tv       = sum(sqrt(U(:).^2 + V(:).^2));
end

function Dx = LinOp(x)
    % Computes horizontal and vertical differences in the geometry
    % component with Dirichlet boundary conditions, corresponds to a linear 
    % operator L. 
    % Inputs -------
    %    x (vector length 2N): vectorized texture and geometry components, each one of size mxn
    % Outputs ------
    %    Dx (matrix size Nx2): first (second) column contains horizontal (vertical) pixel 
    %                              differences in the geometry
    X         = reshape(x(N+1:end),m,n);
    Du        = X;
    Du(:,2:n) = X(:,2:n)-X(:,1:(n-1));   
    Dv        = X;
    Dv(2:m,:) = X(2:m,:)-X(1:(m-1),:);   
    Dx        = [Du(:);Dv(:)];                     
end

function Dtz = LinOpT(z)  
    % Applies the transpose of the linear operator computing horizontal
    % and vertical differences in the geometry image.
    % Inputs -------
    %    z (vector length 2N): vectorized 2D image of size mxn
    % Outputs ------
    %    Dtz (matrix size Nx2): L'*z
    U            = reshape(z(1:N),m,n);                                        
    U(:,1:(n-1)) = U(:,1:(n-1))-U(:,2:n);    
    V            = reshape(z(N+1:end),m,n);                                      
    V(1:(m-1),:) = V(1:(m-1),:)-V(2:m,:);    
    Dtz          = U + V;                                  
    Dtz          = [zeros(N,1);Dtz(:)];
end

function [C_1,C_2,C_3,C_4,C_5,C_6,C_full] = constraints(ut,ug)
    % Computes the inequality constraints.
    % Input -------
    %    ut (vector length N): vectorized texture component
    %    ug (vector length N): vectorized geometry component
    % Output ------
    %    C_1, C_2  (vectors length N): constraints on the minimal and maximal pixel value
    %    C_3, C_4  (vectors length q): constraints on the data fidelity term and noise amplitude
    %    C_5, C_6  (vectors length N): constraints on the texture minimal and maximal value
    %    C_full (vector length 4N+2q): concatenation of all constraints
    Hutug   =  H*(ut+ug);
    C_1     =  ut+ug-xmin;
    C_2     = -ut-ug+xmax;
    C_3     =  Hutug-y(:)+chi;
    C_4     = -Hutug+y(:)+chi;
    C_5     =  ut+xmax/3;
    C_6     =  -ut+xmax/3;
    C_full  =  [C_1;C_2;C_3;C_4;C_5;C_6];
end

end