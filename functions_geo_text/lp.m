
%==========================================================================
% Code provided by the authors of
% Boyd, Stephen, and Lieven Vandenberghe. Convex optimization. 
% Cambridge university press, 2004.
% https://web.stanford.edu/~boyd/cvxbook/cvxbook_examples/chap11/
%
% The code was only modified to stop when a strictly feasible point is 
% found.
%==========================================================================

% Barrier method  for a small LP. Boyd & Vandenberghe, Convex 
% Optimization

function [x] = lp(A, b, c, x0, mu, tol,HH,yy,chi)

% [x, inniters, gaps] = lp(A, b, c, x0, mu, tol)
%
%      minimize    c'*x 
%      subject to  A*x <= b
%
%      maximize    -b'*z
%      subject to  A'*z + c = 0
%                  z >= 0
%
% Barrier method for solving LP within absolute accuracy tol, starting
% with initial t = 1, at a strictly feasible x0.  We assume the 
% problem is strictly dual feasible.
% 
% inniters:  array with number of Newton iters per outer iteration
% gaps:   array with duality gaps at the end of each outer iteration
 

MAXITERS = 100;   
ALPHA = 0.01;
BETA = 0.5;
NTTOL = 1e-5;     % stop inner iteration if lambda^2/2 < NTTOL

[m,n] = size(A);
t = 1;
x = x0;
v=[];
    A1=A*ones(n,1);
    Atildet=(A.*A1)';

for k=1:MAXITERS
    d = b-A*x;
    val = t*c'*x - sum(log(d));
    g = t*c + A'*(1./d);
    D=1./(d.^2);
    %H = A'*diag(1./(d.^2))*A;
    %v = -H\g;  
    H = @(u) A'*(D.*(A*u));
    M=sparse(1:n,1:n,Atildet*D);
    [v,info] = pcg(H,-g,1e-4,2000,M,[],v);  
    fprime = g'*v;
    s = 1;
    while (min(b-A*(x+s*v)) < 0),  s = BETA*s;  end;
    while (t*c'*(x+s*v) - sum(log(b-A*(x+s*v))) > val + ALPHA*s*fprime);s = BETA*s; end;
    x = x+s*v;
    %fprintf('iter %d -fprime/2 %d TOL %d\n',k,-fprime/2,NTTOL)
    if mod(k,10)==0
        fprintf('iteration %d: negative if strictly feasible: %.4e\n', k, max(abs(HH*x(2:end)-yy))-chi)
    end
    if max(abs(HH*x(2:end)-yy))<chi
        fprintf('iteration %d: negative if strictly feasible: %.4e\n', k, max(abs(HH*x(2:end)-yy))-chi)
        return 
    end
    
    if ((-fprime/2) < NTTOL) 
        z = (1./d) .* (1 + (A*v)./d) / t;
        gap = (b-A*x)'*z;  
        %fprintf('iter %d gap = %d\n',k,gap)
        if (gap < tol), return;  end; 

        t = min(t*mu, (m+1)/tol);  
   end;
end;
disp(['Maxiters (', int2str(MAXITERS), ') exceeded.']);