% importData(X,folderName,tag,signaRange,Fs,channelLabels)

% Inputs
% X: A 3-D array of signals of size LxMxN, where L = signal length, M =
% Number of trials, N = number of Channels

% folderName:  The folder where MP directories are created. 
% Default: data folder in the current directory 

% tag: An additional string to distinguish the data folder. 
% Folders are created in foldername/tag. Default: 'test' 

% signalRange: the portion of signal to be analysed. If not specified, the
% entire signal L is analyzed. 

% Fs: Sampling Rate (default 1000)
% channelLabels: The channel numbers of the N channels. Default: 1:N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supratim Ray, 2008 
% Distributed under the General Public License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updates
% 1/10/14: Minor name changes in some programss
% 29/3/16: Minor changes with slashes/backslashes

function importData(X,folderName,tag,signalRange,Fs,channelLabels)

sizeX = size(X);
if length(sizeX) == 2
    % just one Channel;
    L = sizeX(1); N = 1; % M = sizeX(2); 
else
    L = sizeX(1); N = sizeX(3); % M = sizeX(2); 
end

if ~exist('folderName','var')
    folderName = fullfile(pwd,'data');
end

if ~exist('tag','var'),                  tag = 'test';                   end
if ~exist('signalRange','var'),          signalRange = [1 L];            end
if ~exist('Fs','var'),                   Fs = 1000;                      end
if ~exist('channelLabels','var'),        channelLabels = 1:N;            end

makeDirectoryMPP(folderName);

writeMPfiles(X,folderName,tag,signalRange);
[EDF, goodChannels, numTrials] = getEDF(X,Fs,channelLabels);

signalLength = signalRange(2)- signalRange(1)+1;
IDnum=1; firstName=tag;
writeheaderfile(folderName,tag,EDF,goodChannels,IDnum,firstName,numTrials,signalLength);

end