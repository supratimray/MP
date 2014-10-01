% This program changes the file delimiters to slash or backslash depending
% on the platform. No changes are made if the file already has the correct delimiter

function fn = platformSpecificNameMPP(fn)

if ispc
    fn(fn=='/')='\';
else
    fn(fn=='\')='/';
end