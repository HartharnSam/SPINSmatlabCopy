function spins_writer(filename, var)
% Writes the MATLAB matrix var that is indexed as (x,y,z) or (x,z) 
% to a SPINS format file.
% The SPINS format file may be read using spins_reader.

% write to file
fid = fopen(filename, 'wb');
fwrite(fid, permute(var, [2 3 1]), 'double');
fclose(fid);
