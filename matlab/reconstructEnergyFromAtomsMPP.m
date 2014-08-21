% [sumEnergy] = reconstructEnergyFromAtomsMPP(gaborData,N,wrap,atomList)

% Given parameters of a gabor book, this program generates the
% time-frequency map. 

% Inputs
% N - length of the signal
% gaborData(1,:) - atom octave
% gaborData(2,:) - atom frequency (0 to N/2)
% gaborData(3,:) - atom time (0 to N-1)
% gaborData(4,:) - atom modulus
% gaborData(5,:) - atom phase (not used)

% wrap - set to 1 if you want wrapping in time
% atomList - list of atom numbers to be reconstructed. Set to [] if you
% want full reconstruction from all the atoms.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supratim Ray, 2008 
% Distributed under the General Public License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [sumEnergy] = reconstructEnergyFromAtomsMPP(gaborData,N,wrap,atomList)

if ~exist('atomList','var'),         atomList=[];                        end

if isempty(atomList)
    atomList = 1:size(gaborData,2);   
end

limitFreq=1;
sumEnergy = zeros(N/2+1,N);
maxOct = nextpow2(N);

for i=1:length(atomList)
    j=atomList(i);
    oct = gaborData(1,j);
    ksi = gaborData(2,j);
    u   = gaborData(3,j);
    mod = gaborData(4,j);
    
    mod2 = mod^2;
    
    if oct==0                       % Dirac
        sumEnergy(2:N/2,u+1) = mod2*(2/N) + sumEnergy(2:N/2,u+1);
        sumEnergy(1,u+1)     = mod2*(1/N) + sumEnergy(1,u+1);
        sumEnergy(N/2+1,u+1) = mod2*(1/N) + sumEnergy(N/2+1,u+1);
        
    elseif oct==maxOct              % Fourier
        sumEnergy(ksi+1,:) = mod2*(1/N) + sumEnergy(ksi+1,:);
        
    else                            %Gabor
        % Get Time Range
        % Note that the energy decays at exp(-2*pi*((t-u)/s)^2), and hence
        % becomes less than 10^(-11) outside u +- 2s. We exploit this to
        % reconstrcut only part of the signal.
        
        nFFT=N;
        if oct==maxOct-1 % range of +-N
            range = [0 N-1];
            if wrap
                [E0,f] = wignerVille(oct,u,ksi,nFFT,range,limitFreq);
                E1 = wignerVille(oct,u-N,ksi,nFFT,range,limitFreq);
                E2 = wignerVille(oct,u+N,ksi,nFFT,range,limitFreq);
                E = (E0+E1+E2);
                normE = sum(sum(E));
                sumEnergy(f+1,:) = mod2*E/normE + sumEnergy(f+1,:);
            else
                [E,f,normE2] = wignerVille(oct,u,ksi,nFFT,range,limitFreq);
                sumEnergy(f+1,:) = mod2*E/normE2 + sumEnergy(f+1,:);
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
                    [E0,f,normE20] = wignerVille(oct,u,ksi,nFFT,r0,limitFreq);
                    r1 = [N-1-(s2-u) N-1];
                    [E1,~,normE21] = wignerVille(oct,u+N,ksi,nFFT,r1,limitFreq);
                    normE2 = normE20+normE21;
                    
                    sumEnergy(f+1,1+(r0(1):r0(2))) = mod2*E0/normE2 + sumEnergy(f+1,1+(r0(1):r0(2)));
                    sumEnergy(f+1,1+(r1(1):r1(2))) = mod2*E1/normE2 + sumEnergy(f+1,1+(r1(1):r1(2)));
                    
                elseif (u > N-1-s2)
                    r0 = [u-s2 N-1];
                    [E0,f,normE20] = wignerVille(oct,u,ksi,nFFT,r0,limitFreq);
                    r1 = [0 u+s2 - (N-1)];
                    [E1,~,normE21] = wignerVille(oct,u-N,ksi,nFFT,r1,limitFreq);
                    normE2 = normE20+normE21;
                    
                    sumEnergy(f+1,1+(r0(1):r0(2))) = mod2*E0/normE2 + sumEnergy(f+1,1+(r0(1):r0(2)));
                    sumEnergy(f+1,1+(r1(1):r1(2))) = mod2*E1/normE2 + sumEnergy(f+1,1+(r1(1):r1(2)));
                else
                    r = [u-s2 u+s2];
                    [E,f,normE2] = wignerVille(oct,u,ksi,nFFT,r,limitFreq);
                    sumEnergy(f+1,1+(r(1):r(2))) = mod2*E/normE2 + sumEnergy(f+1,1+(r(1):r(2)));
                end
            else
                r = [max(0,u-s2) min(N-1,u+s2)];
                [E,f,normE2] = wignerVille(oct,u,ksi,nFFT,r,limitFreq);
                sumEnergy(f+1,1+(r(1):r(2))) = mod2*E/normE2 + sumEnergy(f+1,1+(r(1):r(2)));
            end
        end
    end
end

end