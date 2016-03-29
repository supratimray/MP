% This program writes the header and the list files
% This program is a gerelization of the writeheaderfile program used to
% write header and list files for the ECoG data. To preserve formatting,
% the same structure is used for any file. Use the function getEDF to get
% the proper EDF and goodchannels variables for your data. IDnum and First_Name
% are optional variables. Numtrials is the number of trials

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supratim Ray, 2008 
% Distributed under the General Public License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function writeheaderfile(foldername,L,EDF,goodchannels, IDnum, First_Name, NumTrials,SignalLength)

if nargin < 7
	NumTrials = input('Enter number of trials: ');
	SignalLength = input('Enter signal length: ');
elseif nargin < 8
	SignalLength = input('Enter signal length: ');
end

labels = EDF.Label;
for i = 1:length(goodchannels)
    tmpL = labels(goodchannels(i),:);
    d = find(tmpL(3:7) == ' ', 1, 'last')+3;
    allL(i,:) = tmpL(d:d+2); %#ok<*AGROW>
    allnum(i,:) = str2num(tmpL(d+3:d+6)); %#ok<*ST2NM>
end

% Open the file
fn = fullfile(foldername,L);
makeDirectoryMPP(fn);

fn1 = fullfile(fn,'ImportData_SIG');
makeDirectoryMPP(fn1);

filename = fullfile(fn1,'sig.hdr');
fp = fopen(filename,'w');

% Taking Data to generate file
Year = '08';
%IDnum = input('Enter the last three digits of the ID: ');
ID = [Year,num2str(IDnum)];
%First_Name = input('Enter the first name: ','s');
Last_Name = ['PY',Year,'N',num2str(IDnum)];

changedefaults = 0;
if changedefaults == 1
    Bytes_per_var = input('Enter Bytes_per_var value: ');
    Name_template = input('Enter Name_template: ');
else
    Bytes_per_var = 8;
    Name_template = 'sig.dat';
end

Samp_rate = EDF.SampleRate(1);
Numb_chans = length(goodchannels);

% Writing the file
fprintf(fp,'%s\n','%PATIENT_INFO');
fprintf(fp,'%s%s\n','ID=',ID);
fprintf(fp,'%s%s\n','First_Name=',First_Name);
fprintf(fp,'%s%s\n','Last_Name=',Last_Name);
fprintf(fp,'%s\n','DOB=01/01/01');

fprintf(fp,'%s\n','%FILE_FMT');
fprintf(fp,'%s%s\n','Samp_rate=',num2str(Samp_rate));
Bytes_per_win = Numb_chans*Bytes_per_var;
fprintf(fp,'%s%d\n','Bytes_per_win=',Bytes_per_win);
fprintf(fp,'%s%s\n','Numb_chans=',num2str(Numb_chans));
fprintf(fp,'%s\n','Byte_order=Little_Endian');
fprintf(fp,'%s\n','File_format=double');
fprintf(fp,'%s%d\n','Bytes_per_var=',Bytes_per_var);
fprintf(fp,'%s\n','Numb_vars=1');
fprintf(fp,'%s\n','Data_offset=0');
fprintf(fp,'%s\n','Chans_per_file=-1');
fprintf(fp,'%s%f\n','Win_len=',(1/Samp_rate));
fprintf(fp,'%s%f\n','Win_shift=',(1/Samp_rate));
fprintf(fp,'%s%s\n','Name_template=',Name_template);
fprintf(fp,'%s\n','Start_chan_no=0');
fprintf(fp,'%s\n','List_file=sig.lst');
fprintf(fp,'%s\n','%CHANNELS');
fprintf(fp,'%s%d\n','Chan_labels=',Numb_chans);

for i = 1:Numb_chans
    fprintf(fp,'%s %d\n',allL(i,:),allnum(i));
end

fprintf(fp,'%s\n','AVG_definitions=0');
fprintf(fp,'%s\n','Chan_units=0');
fprintf(fp,'%s\n','%VARIABLES');
fprintf(fp,'%s\n','Var_labels=0');
fprintf(fp,'%s\n','AVG_definitions=0');
fprintf(fp,'%s\n','Var_units=0');
%fprintf(fp,'%s\n','uV');

calibrate = 0;
if calibrate == 1
    
    fprintf(fp,'%s\n','%CALIBRATION');
    fprintf(fp,'%s%d\n','Chan_calib=',Numb_chans);

    for i=1:Numb_chans
        j = goodchannels(i);
        fprintf(fp,'%1.7f\n',1000/(max(abs(EDF.DigMin(j)),EDF.DigMax(j))));
    end
    fprintf(fp,'%s\n','Var_calib=0');
end

fclose(fp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we write the list file
filename = fullfile(fn1,'sig.lst');
fp = fopen(filename,'w');
fprintf(fp,'%s\t%s\t%f\t%s\n',Name_template, 'Sat 1 Jan 2000 0:00',SignalLength*NumTrials/Samp_rate, '0.0000');
fclose(fp);

end