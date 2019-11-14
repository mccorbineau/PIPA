function wrc=tab2cell(wc,jm,Ml,Mc,N)

%mise en forme des coefs sous la forme w{j}{ml}{mc}
%au départ sont sous la forme [[vert horz diag]' vert horz diag app]'
ind=1;
for j=1:jm,
    for ml=1:Ml,
        for mc=1:Mc,
            %[j,ml,mc]
            %approx
            if j==jm & ml==1 & mc==1
                inda=3*jm+1;
                wrc{j}{ml}{mc}=wc((inda-1)*N+1:inda*N,:);
                %coefs redecomposés
            elseif j~=jm & ml==1 & mc==1
                wrc{j}{ml}{mc}=[];
                %le reste
            else
                wrc{j}{ml}{mc}=wc((ind-1)*N+1:ind*N,:);
                ind=ind+1;
            end
        end
    end
end