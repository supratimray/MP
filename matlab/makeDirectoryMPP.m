% makeDirectoryMPP(foldername)
% Makes the folder if is does not exist.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supratim Ray, 2008 
% Distributed under the General Public License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updates
% 1/10/14: Minor name changes in some programs

function makeDirectoryMPP(fn)

if isdir(fn)==0
    disp(['Creating directory ',fn]);
    mkdir(platformSpecificNameMPP(fn));
end