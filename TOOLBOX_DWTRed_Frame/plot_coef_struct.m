function plot_coef_struct(w)

J=length(w);
Ml=length(w{1});
Mc=length(w{1}{1});

for j=1:J,
    coef=[];
    for ml=1:length(w{j}),
        temp=[];
        for mc=1:length(w{j}{ml}),
            if isempty(w{j}{ml}{mc})
                inter=zeros(size(w{j}{ml+1}{mc+1}));
                temp=[temp,inter];
            else
                temp=[temp,w{j}{ml}{mc}];
            end
        end
        coef=[coef;temp];
    end
    [lc,cc]=size(coef);
    s=1;
    figure
    for k=1:Ml,
        for l=1:Mc,
            tit=[inputname(1),'\{',int2str(j),'\}\{',int2str(k),'\}\{',int2str(l),'\}'];
            subplot (Ml,Mc,s)
            imagesc(coef((k-1)*lc/Ml+1:k*lc/Ml,(l-1)*cc/Mc+1:l*cc/Mc))
            set(gca,'XTick',[1,size(coef((k-1)*lc/Ml+1:k*lc/Ml,(l-1)*cc/Mc+1:l*cc/Mc),2)])
            set(gca,'XTickLabel',{'1',int2str(size(coef((k-1)*lc/Ml+1:k*lc/Ml,(l-1)*cc/Mc+1:l*cc/Mc),2))})
            set(gca,'YTick',[1,size(coef((k-1)*lc/Ml+1:k*lc/Ml,(l-1)*cc/Mc+1:l*cc/Mc),1)])
            set(gca,'YTickLabel',{'1',int2str(size(coef((k-1)*lc/Ml+1:k*lc/Ml,(l-1)*cc/Mc+1:l*cc/Mc),1))})
            axis 'equal' 'tight'
            colormap(gray(256))
            xlabel(tit)
            s=s+1;
        end
    end
end
