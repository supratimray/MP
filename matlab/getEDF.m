% [EDF, NumTrials, goodchannels] = getEDF(X,Fs,channelLabels)

% X is a 3-D array of data as defined in importData.m
% Fs: Sampling frequency
% channelLabels: Is a set of integers to specify the channel numbers.
% Default: 1:number_of_channels. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supratim Ray, 2008 
% Distributed under the General Public License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [EDF, goodChannels, numTrials] = getEDF(X,Fs,channelLabels)

sizeX = size(X);

if length(sizeX) < 3   % only one channel
    numTrials = sizeX(2); numChans = 1; 
else
    numTrials = sizeX(2); numChans = sizeX(3);
end

if nargin < 3
    % Choose default ChannelLabels
    channelLabels = 1:numChans;
else
    % check if the NumChans and length(goodchannels) are the same
    if numChans ~= length(channelLabels)
        error('Length of goodchannels different from the Number of channels specified in X');
    end
end

% We need two parameters from EDF: EDF.SampleRate and EDF.Labels
EDF.SampleRate = Fs;

for i=1:numChans
    EDF.Label(i,:) = repmat(' ',1,16);
    tmpLabel = ['SIG CHN ' num2str(channelLabels(i))];
    EDF.Label(i,1:length(tmpLabel)) = tmpLabel;
end

goodChannels = 1:numChans;

end