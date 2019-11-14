function [Y,S,nRow,nCol,Xtrue,materials] = load_data()
%==========================================================================
% Creates simulated hyperspectral data from the Urban dataset, includes
% artificial atmospheric attenuation.
%
% Output ------------------------------------------------------------------
%    Y               (nB*nPix): hyperspectral data for each pixel and every
%                               wavelength
%    S               (nB*nEnd): library of the spectral signatures of the 
%                               endmembers
%    nRow                (int): number of rows
%    nCol                (int): number of columns
%    Xtrue  (matrix nEnd*nPix): ground-truth of the abundance maps 
%    materials    (cell array): names of the studied materials (endmembers)
%==========================================================================

rng(5); % noise seed, for reproducibility
disp('------------------------------------------------------------------------')
disp('Create noisy data')
disp('------------------------------------------------------------------------')

%%%%%%%% load data
load('data_hyperspectral/Urban')
Xtrue     = Urban.Xtrue;     % ground-truth of the abundance maps
nEnd      = Urban.nEnd;      % number of materials or endmembers
materials = Urban.materials; % list of materials names 
nRow      = Urban.nRow;      % number of rows in the abundance maps
nCol      = Urban.nCol;      % number of columns in the abundance maps
S         = Urban.lib;       % library including the spectral signatures of 
                             % the materials
fprintf('dataset:                  Urban\n')
fprintf('number of materials:      %d\n',nEnd)
fprintf('image size:               %d x %d\n',nRow,nCol)

%%%%%%%% generate attenuated noisy hyperspectral data
%%% add attenuation to satisfy the constraint <=1
%%% attenuation is different for every pixel
eta   =  0.05;          % eta = 0 -> no attenuation
t1    = (1:nRow)/nRow;
t2    = (1:nCol)/nCol;
Latt  =  rand(2,1);
Att   =  kron(exp(-10*(t2-Latt(1)).^2),exp(-10*(t1-Latt(2)).^2)')*5;
Att   =  Att(:)';
Xatt  = (1-eta*Att).*Xtrue; % attenuated abundance maps
Y     =  S*Xatt;            % hyperspectral data

%%%%%%%% add noise
noise_snr = 1.5;                                  % signal-to-noise ratio
noise_std = mean(std(Y,[],2))*10^(-noise_snr/20); % Gaussian noise standard deviation
Y         = Y + randn(size(Y)).*noise_std;
fprintf('noise standard deviation: %.3e\n',noise_std);
end