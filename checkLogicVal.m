function out = checkLogicVal(in)

if isnumeric(in)
	out = (in ~= 0);
elseif ischar(in)
	if strcmpi(in(1), 'y')
		out = true;
	else
		out = false;
	end
else
	error('%s: cannot process value for in: %s', mfilename, class(in));
end
