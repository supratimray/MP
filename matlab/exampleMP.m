% exampleMP: A simple example to illustrate the MP program.
% Supratim Ray, 2008
% just type 'exampleMP' to run this program. 

% Input
% Suppose the input data consists of 3 trials for 2 different conditions,
% as described below.

Fs = 1000; % Sampling frequency
L  = 1024; % signal length
t  = (0:L-1)/Fs;

% condition 1, some gabor blobs in gaussian noise
Y1=gauspuls(t-0.5,30,0.25) + gauspuls(t-0.25,200,0.2) + gauspuls(t-0.75,80,0.5);
inputSignal(:,1,1) = Y1 + 0.01*randn(1,L); % mostly signal
inputSignal(:,2,1) = Y1 + 0.1*randn(1,L);
inputSignal(:,3,1) = Y1 + 1*randn(1,L); % mostly noise

% condition 2, a chirp signal+ a sinusoid
Y2=chirp(t,100,0.5,200) + sin(2*pi*50*t);
inputSignal(:,1,2) = Y2 + 0.01*randn(1,L);
inputSignal(:,2,2) = Y2 + 0.1*randn(1,L);
inputSignal(:,3,2) = Y2 + 1*randn(1,L);


% specify folderName and tag
folderName = [pwd '/data/'];
tag = 'test/';

% Import the data
X = inputSignal;
signalRange = [1 L]; % full range
importData(X,folderName,tag,signalRange,Fs);

% perform Gabor decomposition
Numb_points = L; % length of the signal
Max_iterations = 100; % number of iterations
runGabor(folderName,tag,Numb_points, Max_iterations);


%%%%%%%%%%% Retrieve and display
% Retrieve information
% Note, for the retrieval to work, the folder viewData should be in the
% maltab path. Or else, simply cd to viewData.

disp('cd to viewdata or set viewData in your matlab path');

% plot Trial number 
trialNum=2;

gaborInfo{1} = getGaborData(folderName,tag,1);
gaborInfo{2} = getGaborData(folderName,tag,2);

% reconstruct signal
wrap=1;
atomList=[]; % all atoms

if isempty(atomList)
    disp(['Reconstructing trial ' num2str(trialNum) ', all atoms']);
else
    disp(['Reconstructing trial ' num2str(trialNum) ', atoms ' num2str(atomList(1)) ':' num2str(atomList(end))]);
end

rSignal1 = reconstructSignalFromAtomsMPP(gaborInfo{1}{trialNum}.gaborData,L,wrap,atomList);
rSignal2 = reconstructSignalFromAtomsMPP(gaborInfo{2}{trialNum}.gaborData,L,wrap,atomList);

% reconstruct energy
rEnergy1 = reconstructEnergyFromAtomsMPP(gaborInfo{1}{trialNum}.gaborData,L,wrap,atomList);
rEnergy2 = reconstructEnergyFromAtomsMPP(gaborInfo{2}{trialNum}.gaborData,L,wrap,atomList);
f = 0:Fs/L:Fs/2;

subplot(421);
plot(t,inputSignal(:,trialNum,1)); 
title(['Condition 1, trial: ' num2str(trialNum)]); xlabel('Time (s)');
axis tight

subplot(423);
plot(t,rSignal1,'k');
title('Reconstruction'); xlabel('Time (s)');
axis tight

subplot(223);
pcolor(t,f,rEnergy1); shading interp;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
axis tight

subplot(422);
plot(t,inputSignal(:,trialNum,2));
title(['Condition 2, trial: ' num2str(trialNum)]); xlabel('Time (s)');
axis tight

subplot(424);
plot(t,rSignal2,'k');
title('Reconstruction'); xlabel('Time (s)');
axis tight

subplot(224);
pcolor(t,f,rEnergy2); shading interp;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
axis tight