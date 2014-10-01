% strout = removeIfPresentMPP(strin,tag)
% Removes tag from strin if it is present

function strout = removeIfPresentMPP(strin,tag)

if ~isempty(strin)
    if ~strcmp(strin(end-length(tag)+1:end),tag)
        strout = strin;
    else
        strout = strin(1:end-length(tag));
    end
else
    strout = strin;
end

return