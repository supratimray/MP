% This function computes the Wigner-Ville distribution of an Gabor atom as
% described in WignerVille.doc/pdf in this directory.

% Inputs:

% oct       - log2(scale);               % between 1 & log2(N)-1
% u         - position                   % between 0 & nFFT-1
% ksi       - frequency                  % between 0 & N-1
% nFFT      - size of FFT
% range     - [r1 r2]                    % time range between which WV should be computed
% both r1, r2 should be between 0.. N-1
% Note that the signal length N is not required

% optional variable limitFreq: set to 1 if you want the frequency values
% only in a smaller range in which they are > ~10^(-10). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supratim Ray, 2008 
% Distributed under the General Public License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [E,f,normE2] = wignerVille(oct,u,ksi,nFFT,range,limitFreq)

if ~exist('limitFreq','var'),           limitFreq=1;             end

if ksi > nFFT/2
    disp( 'ksi should not exceed nFFT/2,taking mirror image');
    ksi = nFFT - ksi;
end

s = 2^oct;

% Time values
timeVals = range(1):range(2);
t = (timeVals-u)/s;
timeAxis = exp(-2*pi*t.*t);

if oct==1
    f = 0:nFFT/2;   
    freqAxis = WV(s*(f-ksi)/nFFT) + WV(s*(f+ksi)/nFFT) + WV(s*(f-nFFT+ksi)/nFFT) + ...
        WV(s*(f-nFFT-ksi)/nFFT) + WV(s*(f+nFFT-ksi)/nFFT);

elseif oct==2
    f = 0:nFFT/2;
    freqAxis = WV(s*(f-ksi)/nFFT) + WV(s*(f+ksi)/nFFT) + WV(s*(f-nFFT+ksi)/nFFT);
  
else
    
    if ksi < 2*nFFT/s
        f = 0:nFFT/2;
        freqAxis = WV(s*(f-ksi)/nFFT) + WV(s*(f+ksi)/nFFT);
    elseif ksi > nFFT*(1/2 - 2/s)
        f = 0:nFFT/2;
        freqAxis = WV(s*(f-ksi)/nFFT) + WV(s*(f-nFFT+ksi)/nFFT);
    else
        if limitFreq
            f = ksi-2*nFFT/s:ksi+2*nFFT/s;
        else
            f = 0:nFFT/2;
        end
        freqAxis = WV(s*(f-ksi)/nFFT);
    end
end

E = freqAxis*timeAxis;

% Get normalization constant
if length(f)<nFFT/2+1
    normE2 = sum(sum(E));    % overall energy, square of the norm
else
    % 0 and N/2 frequencies are halved - this is like doubling the remaining
    % frequencies. Normalization takes care of the scales
    E(1,:)=E(1,:)/2; 
    E(end,:) = E(end,:)/2;
    
    normE2 = sum(sum(E));
end

end

function freqAxis = WV(f)
freqAxis = exp(-2*pi*f.*f)'; % equation 60 in mallat and Zhang
end
