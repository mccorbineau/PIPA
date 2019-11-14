disk=fspecial('disk',61);
disk(disk<max(disk(:))/2)=0;
disk(disk>0)=1;
disk_full=zeros(128,128);
disk_full(3:125,3:125)=disk;

length=415;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath Agaricus_bisporus/8bitTIFF

agari=imread('agari0330.tif');
agari_rescale=imresize(agari(420:420+length,420:420+length),[128,128]);
disk_agari = disk_full.*im2double(agari_rescale);

figure;
subplot(131);imagesc(agari);colormap gray
subplot(132);imagesc(agari_rescale);colormap gray
subplot(133);imagesc(disk_agari);colormap gray

%%%%%%%%%%%%%%%%%%%%%%

load('image1024.mat');

glass_rescale=imresize(I(286:286+length,87:87+length),[128,128]);
disk_glass = disk_full.*glass_rescale;

figure;
subplot(131);imagesc(I);colormap gray
subplot(132);imagesc(glass_rescale);colormap gray
subplot(133);imagesc(disk_glass);colormap gray
