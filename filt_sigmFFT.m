% Created on 4 Mar 2015
% Primary author: Kiefer Forseth
%                 kjforseth@gmail.com
%                 602.531.6430
%
% Summary
%   Filter traces into a band with the option of applying an AC notch.
%   NOTE: The AC notch is not presently validated
% 
% Input
%   band: bandpass range (set to [1 300] if you want to just apply notch)
%   notch: logical option, true applies AC filter
%   fs: sampling rate from nkdata
%   fb: filter buffer (trimmed at end; set to 0 if you want no trimming)
%   varargin: add trace matrices of size (# trials) x (# samples)
%
% Output
%   varargout: filtered trace matrices of size (# trials) x (# samples - 2*fb)
% 
% Example
%   Extract Mid-Gamma from single trace matrix:
%       ftraces = filt_sigmFFT([61 119],false,fs,fb,traces;
%   Apply notch only to 2 trace matrices:
%       [ftraces,fbase] = filt_sigmFFT([1 300],true,fs,fb,traces,base);
%
% Change log
%   7 July 15: fixed AC notch

function [varargout] = filt_sigmFFT(band,notch,fs,fb,varargin)

% handle inputs
if iscell(band)
    rolloff = band{2};
    band = band{1};
else
    rolloff = [3 3];
end

% sigmoid inline functions
sig.h = @(z, cutoff, rolloff) 1 ./ (1 + exp(rolloff*(+cutoff - z)));
sig.l = @(z, cutoff, rolloff) 1 ./ (1 + exp(rolloff*(-cutoff + z)));

% filter each input signal
varargout = cell(1,length(varargin));
for n = 1:length(varargin)
    % filter parameters
    nfft = 2^nextpow2(size(varargin{n},2));
    f = fs/2*linspace(0,1,nfft/2);
    
    % design filter
    d = [sig.l(f,band(2),rolloff(2)).*sig.h(f,band(1),rolloff(1)) zeros(1,nfft/2)];
    d = 2*d/max(d);
    
    if notch
        ACbands = bsxfun(@plus,60*(1:8)'*[1 1],[-0.5 +0.5]);
        
        AC = zeros(1,nfft);
        for b = 1:size(ACbands,1)
            side = sig.l(f,ACbands(b,1),5).*sig.h(f,ACbands(b,2),5);
            side = (side/max(side)); % normalize here in case bands of unequal width
            AC = max([AC; [side fliplr(side)]]);
        end
        AC = 1 - (AC/max(AC));
        
        d = d.*AC; % bandpass & notch filters
    end
    
    % apply filter
    F = fft(varargin{n},nfft,2);
    F = bsxfun(@times,F,d);
    tmp = ifft(F,[],2);
    varargout{n} = tmp(:,fb+1:size(varargin{n},2)-fb);
end