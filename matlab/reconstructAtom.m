% [signal,t,normSignal] = reconstructAtom(oct,u,ksi,fi,N,range)

% This program reconstructs the signal given the properties of the gabor atom

% Inputs:
% N - length of signal
% oct       - log2(scale);               % between 1 & log2(N)-1
% u         - position                   % between 0 & N-1
% ksi       - frequency                  % between 0 & N-1
% fi        - angle                      % in radians
% range - [r1 r2]                        % both r1, r2 between 0.. N-1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supratim Ray, 2008 
% Distributed under the General Public License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [signal,t,normSignal] = reconstructAtom(oct,u,ksi,fi,N,range)

if ~exist('range','var'),          range=[];           end

if isempty(range)
    t=0:N-1;
else
    if ((range(1) < 0) || range(2) > N-1)
        error('limits out of range. Change the range variable');
    else
        t=range(1):range(2);
    end
end

s = 2^oct;
cosTerm = cos(2*pi*ksi*t/N+fi);
expTerm = G((t-u)/s);
signal = cosTerm.*expTerm;
normSignal = sqrt(sum(signal.*signal));
end

function S = G(t)
    S = exp(-pi*t.*t);  % The constant 2^(1/4) is not needed because normalization takes care of it
end