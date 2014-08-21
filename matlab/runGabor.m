% runGabor(foldername,tag,Numb_points, Max_iterations)

% This program run the MP algorithm.
% foldername:  The folder where MP directories are created. 
% tag: An additional string to distinguish the data folder.  
% Numb_points is the length of the signal in points (default = 1024);
% Max_iterations = maximum number of Gabor atoms for each signal. (default 500)

% Note:
% Please change the sourcefolder variable on line 21 to the folder where
% the executable gabord is located.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supratim Ray, 2008 
% Distributed under the General Public License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function runGabor(foldername,tag,Numb_points, Max_iterations)

sourcefolder = '/Users/sray/Research/timeFrequency/mpp/MP/source/';

if ~exist('Numb_points')        Numb_points=1024;       end
if ~exist('Max_iterations')     Max_iterations=500;     end

foldername=appendIfNotPresent(foldername,'/');
tag=appendIfNotPresent(tag,'/');

fnin = [foldername tag];
fn = [fnin 'ImportData_SIG/GaborMP/'];
makeDirectory(fn);
fnout = [fnin 'GaborMP/'];
makeDirectory(fnout);

Shift_points = Numb_points;

% Get information about the number of channels from the header file
filename = [fnin 'ImportData_SIG/sig.hdr'];
Numb_chans = getFromFile(filename,'Numb_chans');
All_chans = Numb_chans;


% Create the control file (called local.ctl)
filename = [fn 'local.ctl'];
fp = fopen(filename,'w');

fprintf(fp,'%s\n','# Template for running Gabor MP Analysis');
fprintf(fp,'%s\n','%INPUT_OUTPUT');
fprintf(fp,'%s\n','Numb_inputs=1');
fprintf(fp,'%s\n','Numb_outputs=1');
fprintf(fp,'%s\n','Mode=parallel');
fprintf(fp,'%s\n','%INPUT');
fprintf(fp,'%s%s\n','Path=',[fnin 'ImportData_SIG/']);
fprintf(fp,'%s\n','Header_file=sig.hdr');
fprintf(fp,'%s\n','Calibrate=1');
fprintf(fp,'%s%d\n','Numb_points=',Numb_points);
fprintf(fp,'%s%d\n','Shift_points=',Shift_points);
fprintf(fp,'%s\n','%OUTPUT');
fprintf(fp,'%s%s\n','Path=',fnout);
fprintf(fp,'%s%d\n','All_chans=',All_chans);
fprintf(fp,'%s%d\n','Numb_chans=',Numb_chans);
fprintf(fp,'%s\n','Start_chan=1');
fprintf(fp,'%s\n','Start_chan_no=0');
fprintf(fp,'%s\n','Header_file=book.hdr');
fprintf(fp,'%s\n','Type=book');
fprintf(fp,'%s\n','File_format=double');
fprintf(fp,'%s\n','Name_template=mp#.bok');
fprintf(fp,'%s\n','Max_len=600');
fprintf(fp,'%s\n','Chans_per_file=-1');
fprintf(fp,'%s\n','%GABOR_DECOMPOSITION');
%fprintf(fp,'%s\n','#Energy=50& int & &&Table=Gabor_MP&');
%fprintf(fp,'%s\n','#Energy=50& int & &&Table=Gabor_MP&');
fprintf(fp,'%s%d\n','Max_Iterations=',Max_iterations);
%fprintf(fp,'%s\n','#Coherence=1& int & &&Table=Gabor_MP&');
fclose(fp);

% Run the GaborMP program (Call the gabord function)
commandline = [sourcefolder 'gabord ' filename];
unix(commandline)