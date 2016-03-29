% writeMPfiles(X,foldername,tag,range)
% This program is called by importData.m

% Inputs
% X: a 3-D array with dimensions LxMxN, where L = length
% of the signal, M = Number of trials and N = Number of channels

% foldername and tag: the MP data files are created in
% foldername/tag/ImportData/

% range (= [L1 L2]) is an optional two vector array that specifies the chunk of
% signal to be taken for gabor analysis. X(:,L1:L2) is taken for analysis.
% If not specified, the entire signal is taken. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supratim Ray, 2008 
% Distributed under the General Public License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 29/3/16: Minor changes with slashes/backslashes

function writeMPfiles(X,foldername,tag,range)

sizeX = size(X);
if length(sizeX) == 2
    % just one Channel;
    L = sizeX(1); M = sizeX(2); N = 1;
else
    L = sizeX(1); M = sizeX(2); N = sizeX(3);
end

if nargin < 4
    % range not specified. Taking the entire signal
    L1 = 1; L2 = L;
    disp(['Taking default signal length of ' num2str(L)]);
else
    L1 = range(1); L2 = range(2);
    disp(['Taking signal from ' num2str(L1)  ' to ' num2str(L2)]); 
end

fn = fullfile(foldername,tag);
makeDirectoryMPP(fn);

fn1 = fullfile(fn,'ImportData_SIG');
makeDirectoryMPP(fn1);

if N > 1    % More than one channel
    for cNum = 1:N
        if cNum <= 10
            fidw = fopen(fullfile(fn1,['sig.dat.00' num2str(cNum-1)]),'w');
        elseif cNum <= 100
            fidw = fopen(fullfile(fn1,['sig.dat.0' num2str(cNum-1)]),'w');
        else
            fidw = fopen(fullfile(fn1,['sig.dat.' num2str(cNum-1)]),'w');
        end

        for k = 1:M
            TW = X(L1:L2,k,cNum);
            fwrite(fidw,TW,'double');
            clear TW
        end
        fclose(fidw);
    end

else	% Only one channel
    fidw = fopen(fullfile(fn1,'sig.dat.000'),'w');
    if M > 1	% More than 1 trial
        for k = 1:M
            TW = X(L1:L2,k);
            fwrite(fidw,TW,'double');
            clear TW
        end
        fclose(fidw);
    else
        TW = X(L1:L2);
        fwrite(fidw,TW,'double');
        fclose(fidw);
    end
end