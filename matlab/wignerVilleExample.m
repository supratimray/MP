function wignerVilleExample(N,ksi,oct)

centerFreq = ksi/N;
freqVals = -2*N:2*N;
scale = 2^oct;

% Note that for omega = 2*pi*k/N, omega/2*pi is simply k/N, k = 0, ... N-1;
% It is then scaled by scale^2.

WV{1} = exp(-2*pi*scale^2*((freqVals)/N - centerFreq).^2);
WV{2} = exp(-2*pi*scale^2*((freqVals-N)/N - centerFreq).^2);
WV{3} = exp(-2*pi*scale^2*((freqVals+N)/N - centerFreq).^2);

%%%%%%%%%% Mirror image

WV{4} = exp(-2*pi*scale^2*((freqVals)/N + centerFreq).^2);
WV{5} = exp(-2*pi*scale^2*((freqVals-N)/N + centerFreq).^2);
WV{6} = exp(-2*pi*scale^2*((freqVals+N)/N +centerFreq).^2);

sumWV =zeros(1,length(freqVals));

colornames = 'bcgrmy';
    
for i=1:6
    plot(freqVals,WV{i},colornames(i)); hold on
    sumWV =sumWV+WV{i}/2;
end
plot(freqVals,sumWV,'k');

for i=1:6
    plot(freqVals,WV{i},[colornames(i) 'o']);
end
plot(freqVals,sumWV,'ko');

axis tight;
xlabel('Frequency (n)');

legend('w-c','w-N-c','w+N-c','w+c','w-N+c','w+N+c','sum');
axis([-N/2 N 0 1.5])