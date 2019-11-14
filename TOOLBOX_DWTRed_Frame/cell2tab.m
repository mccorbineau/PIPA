function wr=cell2tab(wrc,jm,Ml,Mc)

wr=[];
for j=1:jm,
    for ml=1:Ml,
        for mc=1:Mc,
            if ml~=1 | mc~=1
                wr=[wr;wrc{j}{ml}{mc}];
            end
        end
    end
end
wr=[wr;wrc{jm}{1}{1}];