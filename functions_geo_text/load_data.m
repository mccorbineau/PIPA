function [x_true,y,chi,H,weightTV] = load_data(sample_name)
%==========================================================================
% Creates a noisy sinogram by applying the radon transform on a real image 
% obtained from a high-resolution X-ray CT facility.
%
% Input -------------------------------------------------------------------
%    sample_name (str): name of the sample ('agaricus' or 'glass')
%
% Output ------------------------------------------------------------------
%    x_true (matrix size mxn): ground-truth (used only to compute the 
%                              signal-to-noise ratio)
%    y      (vector length q): observed data
%    chi             (double): measurement uncertainty
%    H      (matrix size qxN): observation operator (Radon transform)
%    weightTV        (double): regularization parameter 
%==========================================================================

rng(3)             % noise seed for reproducibility
disp('------------------------------------------------------------------------')
disp('Create a noisy sinogram')
disp('------------------------------------------------------------------------')

xmax = 1;
load('data_geo_text/H','H') % radon operator

% disk mask to simulate a computed tomography acquisition
disk=fspecial('disk',61);
disk(disk<max(disk(:))/2)=0;
disk(disk>0)=1;
disk_full=zeros(128,128);
disk_full(3:125,3:125)=disk;

length = 415;
switch sample_name
    case 'glass'
        load('data_geo_text/image1024.mat','I');
        weightTV      = 0.25;  % regularization parameter set by hand
        glass_rescale = imresize(I(286:286+length,87:87+length),[128,128]);
        x_true        = disk_full.*glass_rescale;
    case 'agaricus'
        weightTV      = 0.5;   % regularization parameter set by hand
        agari         = imread('data_geo_text/Agaricus_bisporus/8bitTIFF/agari0330.tif');
        agari_rescale = imresize(agari(420:420+length,420:420+length),[128,128]);
        x_true        = disk_full.*im2double(agari_rescale);   
end
%normalization
x_true = xmax.*x_true./max(x_true(:));
% apply observation and degradation model for Computed Tomography
y   = H*x_true(:);                      % tomographic data
chi = 0.02*max(abs(y(:)));              % amplitude of uniform noise 2% of max value
y   = y+chi.*(-1+2.*rand(size(H,1),1)); % noisy sinogram

fprintf('sample:           %s\n',sample_name)
fprintf('image size:       128 x 128\n')
fprintf('noise amplitude:  %.2f\n',chi)
end