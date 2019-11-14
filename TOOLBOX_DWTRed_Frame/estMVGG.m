function [beta,alpha] = estMVGG(z,betac)
%[beta,alpha] = estMVGG(z,betac)

if nargin == 1
   prec = 1e-2;
   betac = prec:prec:2;
end


MV = [];
for beta = betac;
   mom = mean(abs(z).^beta);
   MV = [MV 1/beta-log(beta/gamma(1/beta))+1/beta*log(beta*mom)];
end

[MVmin,imin] = min(MV);

beta = betac(imin);
alpha = (beta*mean(abs(z).^beta))^(1/beta);
