function [snr_global,snr_mats] = my_snr(X,Xtrue)
%==========================================================================
% Computes the global signal-to-noise ratio (SNR) for the abundance maps, 
% along with the SNR for each abundance map (each material).
%
% Input -------------------------------------------------------------------
%    X     (matrix nEnd*nPix): estimated abundance map
%    Xtrue (matrix nEnd*nPix): ground-truth
%
% Output ------------------------------------------------------------------
%    snr_global    (double): global SNR in dB
%    snr_mats (vector nEnd): SNR for each material in dB
%==========================================================================

%%% global signal-to-noise ratio
se_global  = sum((X(:)-Xtrue(:)).^2);
ref_global = sum(Xtrue(:).^2);
snr_global = 10*log10(ref_global./se_global);

%%% SNR for each material
se_mats  = sum((X'-Xtrue').^2);
ref_mats = sum(Xtrue'.^2);
snr_mats = 10*log10(ref_mats./se_mats)';
end