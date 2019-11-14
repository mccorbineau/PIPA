function [lambda1,lambda2] = EstimateParameter(I,lambda,M,rm)

[xcc,H1,H2]=Mdwt2_freq(I,M,rm);

betarand = 1;
Ml=size(H1,1);
Mc=size(H2,1);
temp=[];

for j=1:rm
    for ml=1:Ml,
        for mc=1:Mc,
            tempsb=xcc{j}{ml}{mc};
            %Approx
            if j==rm && ml==1 && mc==1
                mu = mean2(xcc{rm}{1}{1});
                tempsb=xcc{rm}{1}{1}(:)- mu ;
                tempsb=tempsb(:)';
                [p,alphap] = estMVGG(tempsb,betarand);
                lambda1=lambda/alphap^p;
            else
                %Details
                temp=[temp,tempsb(:)'];
            end
        end
    end
end
%Details : computation of lambda : one for all the subbands
[p,alphap] = estMVGG(temp,betarand);
lambda2=lambda/alphap^p;

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