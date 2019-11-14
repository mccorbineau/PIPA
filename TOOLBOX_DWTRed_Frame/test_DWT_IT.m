clear all
close all

%taille de l'image
N=256;
%nb de niveaux max que l'on peut faire
J = log2(N);
%niveau de résolution
jm=2;
%filtres
qmf = MakeONFilter('Symmlet',8);
%cas dyadique
Ml=2;
Mc=2;
%image
%x=phantom(N);
%randn('seed',0)
x=randn(256);

figure
imagesc(x)
colormap(gray)
axis image

%transformée en ondelettes invariance par translation
wc = FWT2_TI(x,J-jm,qmf);

%reconstruction
y = IWT2_TI(wc,J-jm,qmf);

disp('Tight frame ?:');
mu=(norm(wc,'fro').^2/norm(x,'fro').^2)


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

%Est-ce que la norme est conservée ?
%Norme des coefs réarrangés
Nwc = 0;
for r = 1:jm
    for m1 = 1:Ml
        for m2 = 1:Mc
            if ~isempty(wrc{r}{m1}{m2})
                Nwc = Nwc+norm(wrc{r}{m1}{m2},'fro')^2;
            end
        end
    end
end

%différence avec la norme des coefs de depart
%disp('Diff entre coef et coef mise en forme {}{}{} : ')
%norm(wc,'fro')^2-Nwc

%remise sous la forme [[vert horz diag]' vert horz diag app]'
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

%Ici wr et wc doivent être exactement les mêmes
if norm(wc-wr)~=0
    error('Pb dans la mise en forme inverse des coefs');
end

figure
imagesc(wc)
colormap(gray)
axis image

%figure
plot_coef_struct(wrc)
colormap(gray)
axis image

figure
imagesc(y)
colormap(gray)
axis image

figure
imagesc(x-y)
colormap(gray)
colorbar

disp('Norm image dep : ');
norm(x,'fro')^2
disp('Norm image rec : ');
norm(y,'fro')^2
disp('Erreur : ');
norm(x-y,'fro')^2

%ADJOINT
%il faut u de la taille coef
u=[];
for i=1:3*jm+1
    u=[u;randn(size(x))];
end
tiu=IWT2_TI_adj(u,J-jm,qmf);
%il faut v de la taille image
v=randn(size(x));
tv=FWT2_TI(v,J-jm,qmf);

disp('Premier prod scalaire : ')
(u(:)'*tv(:))
disp('Deuxième prod scalaire : ')
(tiu(:)'*v(:))