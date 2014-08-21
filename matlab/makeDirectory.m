% makeDirectory(foldername,isWin)
% Makes the folder if is does not exist.
% isWin: set to 1 if on a windows environment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supratim Ray, 2008 
% Distributed under the General Public License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function makeDirectory(fn,isWin)

if ~exist('isWin','var')      isWin=0;        end

if isWin
    % Convert to slashes if fn has backslashes
    fn(find(fn=='/'))='\';
end

if isdir(fn)==0
    disp(['Creating directory ',fn]);
    if ispc
        % Have to do this separately as mkdir as a dos command cannot have
        % trailing slashes. BK.
        mkdir(fn)
    else
        unixcom = ['mkdir ', fn];
        unix(unixcom);
    end
end