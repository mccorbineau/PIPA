function X = Wav_inv_mult(Wav_inv,WX)
%==========================================================================
% This function applies the inverse wavelet transforms to every channel of 
% WX.
%
% Input -------------------------------------------------------------------
%    WX (matrix P*N): wavelet coefficients of an image with multiple 
%                           channels
%    Wav_inv (function): inverse wavelet transform operator
%
% Output ------------------------------------------------------------------
%    X  (matrix P*N): variable whose wavelet coefficients are WX
%==========================================================================

P = size(WX,1); % number of channels
X = zeros(size(WX));

for p = 1:P 
    X(p,:) = Wav_inv(WX(p,:)); % inverse wavelet transform for each channel
end
end
