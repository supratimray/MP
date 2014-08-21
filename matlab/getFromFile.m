function strval = getFromFile(filename,str)
L = length(str);
X = textread(filename,'%s'); 
strNum = X{find(strncmp(X,str,L)==1)};
strval = str2num(strNum(L+2:end));
end