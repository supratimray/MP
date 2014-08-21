% [signal] = reconstructSignalFromAtomsMPP(gaborData,N,wrap,atomList)
% Given parameters of a gabor book, this program generates the signal

% Inputs

% N - length of the signal
% gaborData(1,:) - atom octave
% gaborData(2,:) - atom frequency (between 0 and N/2)
% gaborData(3,:) - atom time      (between 0 and N-1)
% gaborData(4,:) - atom modulus
% gaborData(5,:) - atom phase

% wrap - set to 1 if you want wrapping (periodization) in time
% atomList - list of atom numbers to be reconstructed. Set to [] if you
% want full reconstruction from all the atoms.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supratim Ray, 2008 
% Distributed under the General Public License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [signal] = reconstructSignalFromAtomsMPP(gaborData,N,wrap,atomList)

if ~exist('atomList','var'),         atomList=[];                        end

if isempty(atomList)
    atomList = 1:size(gaborData,2);  
end

signal = zeros(1,N);
maxOct = nextpow2(N);

for i=1:length(atomList)
    j=atomList(i);
    oct = gaborData(1,j);
    ksi = gaborData(2,j);
    u   = gaborData(3,j);
    mod = gaborData(4,j);
    fi  = gaborData(5,j);
    
    if oct==0                       % Dirac
        signal(u+1) = mod*sign(cos(fi)) + signal(u+1);
        
    elseif oct==maxOct              % Fourier
        
        S = cos(2*pi*ksi*(0:N-1)/N+fi);
        normS = sqrt(sum(S.^2));
        signal = mod*S/normS + signal;
        
    else                            %Gabor
        % Get Time Range
        % Note that the energy decays at exp(-2*pi*((t-u)/s)^2), and hence
        % becomes less than 10^(-11) outside u +- 2s. We exploit this to
        % reconstruct only part of the signal.

        if oct==maxOct-1 % range of +-N
            r = [0 N-1];
            if wrap
                S0 = reconstructAtom(oct,u,ksi,fi,N,r);
                S1 = reconstructAtom(oct,u-N,ksi,fi,N,r);
                S2 = reconstructAtom(oct,u+N,ksi,fi,N,r);
                S = (S0+S1+S2);
                normS = sqrt(sum(S.^2));
                signal = mod*S/normS + signal;
            else
                [S,~,normS] = reconstructAtom(oct,u,ksi,fi,N,r);
                signal = mod*S/normS + signal;
            end
        else
            % range on each side of u is at most 2s < N/2, and hence at
            % most only one side needs to be wrapped. Further, the segments
            % are disjointed, so that the norm can be computed by simply
            % adding the norms of the two segments separately.
            
            s2 = 2^(oct+1);
            
            if wrap
                if (u < s2)
                    r0 = [0 u+s2];
                    [S0,t0,normS0] = reconstructAtom(oct,u,ksi,fi,N,r0);
                    r1 = [N-1-(s2-u) N-1];
                    [S1,t1,normS1] = reconstructAtom(oct,u+N,ksi,fi,N,r1);
                    normS = sqrt(normS0^2+normS1^2);
                    
                    signal(1+t0) = mod*S0/normS + signal(1+t0);
                    signal(1+t1) = mod*S1/normS + signal(1+t1);
                    
                elseif (u > N-1-s2)
                    r0 = [u-s2 N-1];
                    [S0,t0,normS0] = reconstructAtom(oct,u,ksi,fi,N,r0);
                    r1 = [0 u+s2-(N-1)];
                    [S1,t1,normS1] = reconstructAtom(oct,u-N,ksi,fi,N,r1);
                    normS = sqrt(normS0^2+normS1^2);
                    
                    signal(1+t0) = mod*S0/normS + signal(1+t0);
                    signal(1+t1) = mod*S1/normS + signal(1+t1);
                else
                    r = [u-s2 u+s2];
                    [S,t,normS] = reconstructAtom(oct,u,ksi,fi,N,r);
                    signal(1+t) = mod*S/normS + signal(1+t);
                end
            else
                r = [max(0,u-s2) min(N-1,u+s2)];
                [S,t,normS] = reconstructAtom(oct,u,ksi,fi,N,r);
                signal(1+t) = mod*S/normS + signal(1+t);
            end
        end
    end
end

end