% strout = appendIfNotPresent(strin,tag)
% Appends tag to strin unless it is present already

% 1/10/14
% 'MPP' appended to the program name so that it does not overload another
% program by the same name in the CommonPrograms folder

function strout = appendIfNotPresentMPP(strin,tag)

if ~isempty(strin)
    if strcmp(strin(end-length(tag)+1:end),tag)
        strout = strin;
    else
        strout = [strin tag];
    end
else
    strout = strin;
end

return