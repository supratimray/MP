% makeDirectory(foldername,isWin)
% Makes the folder if is does not exist.
% isWin: set to 1 if on a windows environment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supratim Ray, 2008 
% Distributed under the General Public License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function makeDirectory(fn)

if isdir(fn)==0
    disp(['Creating directory ',fn]);
    mkdir(platformSpecificName(fn));
end