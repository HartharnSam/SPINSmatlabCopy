function params = spins_params()
% SPINS_PARAMS   Parses the spins.conf file into a structure.
%
%    parms = spins_params()
%
%    David Deepwell, 2015.

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
    disp('spins.conf was not found.')
    params = struct();
end

% check if grid is mapped
if isfield(params,'mapped_grid') == false
    params.mapped_grid = 'false';
    fprintf('The grid was not declared to be mapped or unmapped.\nIt is now assumed to be unmapped, This WILL make problems if incorrect.\n')
end
