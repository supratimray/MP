% Converts a number to string and appends the necessary zeros
function str = conv2Str(n)

if n<10
	str = ['00' num2str(n)];
elseif n < 100
	str = ['0' num2str(n)];
else 
	str = num2str(n);
end