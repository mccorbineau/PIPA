function WX = Wav_mult(Wav,X)
%==========================================================================
% Computes the wavelet coefficients for every channel of the image X.
%
% Input -------------------------------------------------------------------
%    X  (matrix P*N): variable with multiple channels or endmembers
%    Wav  (function): wavelet transform operator
%
% Output ------------------------------------------------------------------
%    WX (matrix P*N): wavelet coefficients of X
%==========================================================================

P  = size(X,1); % number of endmembers
WX = zeros(size(X)); 

for p = 1:P 
    WX(p,:) = Wav(X(p,:)); % wavelet transform for each channel
end
end
