function Z = conv2Log(Y,Thres)

if nargin==1
    Thres = 10^(-16);
end

%%%%%%%%%%%%%% Replace anything less than Thres with Thres %%%%%%%%%%%%%%%%
Y(Y<Thres) = Thres;
Z=log10(Y);

return