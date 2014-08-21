% strout = appendIfNotPresent(strin,tag)
% Appends tag to strin unless it is present already

function strout = appendIfNotPresent(strin,tag)

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