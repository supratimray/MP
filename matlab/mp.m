% This program was written by Dr. Bart Krekelberg while developing the code
% to be windows compatible. Not used/maintained anymore.

function [energy,frequency,signal,time,gaborInfo] = mp(S,varargin)
% Run matching pursuit spectral analyis on the signal in S
% 
% INPUT
% S = Signal with [nrSamples nrTrials nrChannels] 
%       NOTE: nrSamples should be a power of 2.
%                   
% 'maxIterations' Max iterations for the MP algorithm [100]
% 'samplingFreq' Sampling frequency of the data in S. [1000] Hz
% 'workingDir'  Which directory to use as a working directory [tempdir]
% 'graphics'    Show graphical results or not. [false]
% 
% OUTPUT
% energy    = The power spectrum. [nrFrequencies nrSamples nrTrials nrChannels]
% frequency = Frequency axis [nrFrequencies] in Hz.
% signal    =  The reconstructed signal. [nrSamples nrTrials nrChannels]
% time      =  Time base [nrSamples] in s.
% gaborInfo = Raw output of the MP algorithm. {nrChannels} 
%
% NOTE
% 
% NOTE
% This is a wrapper around the Matlab code that Supratim Ray (sray@cns.iisc.ernet.in) wrote
% with some  minor changes to make it compatible with the Windows executable that I
% compiled.
% 
% 
%
% BK - May 2013.
nout=nargout;
p =inputParser;
p.addParamValue('maxIterations',500,@isnumeric); % 
p.addParamValue('samplingFreq',1000,@isnumeric); % 
p.addParamValue('workingDir',tempdir,@ischar);   %
p.addParamValue('graphics',false,@islogical);    % 
p.parse(varargin{:});


if ischar(S) && isdir(S)
    % Special case to test the Windows code against data that SR analyzed
    % on Linux =- the S is the directory name
    [pth,subDir]= fileparts(S);
    [workingDir,dirName] =fileparts(pth);
    nrSamples= 1024;nrChannels = 27;nrTrials=1;
    pushd;
    cd(workingDir);
else
    workingDir = p.Results.workingDir;
    dirName = './MP'; % All in and output will be saved to files in this directory and its subs.
    subDir    = 'data';  % Data will be saved here temporarily
    [nrSamples,nrTrials,nrChannels] = size(S); 
    if log2(nrSamples) ~=round(log2(nrSamples))
        error(['nrSamples should be a power of 2. It is ' num2str(nrSamples)]);
    end
    % Save the data to files.    
    pushd;
    cd(workingDir);
    importData(S,dirName,subDir,[1 nrSamples],p.Results.samplingFreq);    
end
% Have to use a relative path to store files, without  c:\ in it so that the UNIX based
% C-code does not get confused..but we dont know whether we have write
% permission locally. Solution is to cd into the tempdir.

% Do the Gabor decomposition
runGabor(dirName,subDir,nrSamples, p.Results.maxIterations);
%Retrieve the information.
gaborInfo =cell(1,nrChannels);
for i=1:nrChannels 
    gaborInfo{i} = getGaborData(dirName,subDir,i);
end
% Remove temp files and return to the working directory.
popd;
rmdir(fullfile(workingDir,dirName),'s'); 

% Reconstruct
frequency = 0:p.Results.samplingFreq/nrSamples:p.Results.samplingFreq/2;
time      = (0:nrSamples-1)/p.Results.samplingFreq;
if nout >2
    signal  =NaN(nrSamples,nrTrials,nrChannels);
end
energy =NaN(length(frequency),nrSamples,nrTrials,nrChannels);
wrap=1;
atomList=[]; % all atoms
for c=1:nrChannels
    for t=1:nrTrials        
        if nout>2
            signal(:,t,c) = reconstructSignalFromAtomsMPP(gaborInfo{c}{t}.gaborData,nrSamples,wrap,atomList);
        end
        energy(:,:,t,c) = reconstructEnergyFromAtomsMPP(gaborInfo{c}{t}.gaborData,nrSamples,wrap,atomList);    
    end
end

%%  Graphical output : one figure per channel. All trials in one figure.
% The spectrum is overlaid on top of the signal itself.
if (p.Results.graphics)   
    for c=1:nrChannels
        figure(c)
        clf;
        cntr=0;
        for t=1:nrTrials
            cntr= cntr+1;
            subplot(ceil(sqrt(nrTrials)),floor(sqrt(nrTrials)),cntr);
            if nout>2
                hold on
                h=plot(time*1000,min(frequency)+(max(frequency)-min(frequency)).*(signal(:,t,c)-min(signal(:,t,c)))./(max(signal(:,t,c))-min(signal(:,t,c))),'Color',[0.8 0.8 0.8]);               
            end
            h=imagesc(time*1000,frequency,energy(:,:,t,c));
            set(h,'AlphaData',0.8);
            axis xy;
            axis tight
            colormap hot
            xlabel 'Time (ms)'
            ylabel 'Frequency (Hz)'
            title (['Trial ' num2str(t)]);
        end
        suptitle(['Channel ' num2str(c)])
    end
end
