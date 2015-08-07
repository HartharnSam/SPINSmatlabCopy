function params = spins_params()
%  SPINS_PARAMS   Parses the spins.conf file into a structure.
%
%  Usage:
%    parms = spins_params()
%
%  Inputs:
%    n/a
%
%  Outputs:
%    params	- a structure containing fields and values from spins.conf
%
%  David Deepwell, 2015

try
    fileID = fopen('spins.conf');
    tline = fgetl(fileID);
    while ischar(tline)
        if ~isempty(tline)
            lin = strtrim(strsplit(tline,'='));
            nam = lin{1};
            val = lin{2};
            if is_str_numeric(val);
                val = str2double(val);
            end
            params.(nam) = val;
        end
        tline = fgetl(fileID);
    end
    fclose(fileID);
catch
    warning('spins.conf was not found.')
    params = struct();
end

