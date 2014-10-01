% 1/10/14
% 'MPP' appended to the program name so that it does not overload another
% program by the same name in the CommonPrograms folder

function Z = conv2LogMPP(Y,Thres)

if nargin==1
    Thres = 10^(-16);
end

%%%%%%%%%%%%%% Replace anything less than Thres with Thres %%%%%%%%%%%%%%%%
Y(Y<Thres) = Thres;
Z=log10(Y);

return