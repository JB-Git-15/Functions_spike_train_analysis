function out = bload(fname, arg1, startpos, datasize, skip)
% loads a binary file
%
% usage:  bload(fname, [NumChannels SamplesToLoad], startpos, datasize)
% starts reading from startpos BYTES into the file (not samples)
% datasize defaults to 'short'

if (nargin < 4)
	datasize = 'short';
end

fp = fopen(fname, 'r');

if (nargin >=3)
	  status = fseek(fp, startpos, 'bof');
	 % status = fopen(fp, startpos, 'bof');

    if (status ~= 0)
		error('Error doing fseek');
	end
end

if (nargin >= 5)
	out = fread(fp, arg1, datasize, skip);
else
	out = fread(fp, arg1, datasize);
end;

fclose(fp);