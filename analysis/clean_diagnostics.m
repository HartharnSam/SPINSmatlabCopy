function [new_diagnos] = clean_diagnostics()
% clean_diagnostics reads in the diagnostic file and removed repeated
% times, and saves clean version as a .mat file

files_to_clean = {'diagnostics','stresses_top','stresses_bottom', 'plot_times'}; % without file extensions (all are .txt)
%%
for nn = 1:length(files_to_clean)
    % read analysis file (if it exists
    diag_file = [files_to_clean{nn},'.txt'];
    if exist(diag_file, 'file') == 2
        try
            opts = detectImportOptions(diag_file);

            % Specify range and delimiter
            opts.DataLines = [2, Inf];
            opts.Delimiter = ",";

            % Specify file level properties
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule = "read";
            %opts.ReadVariableNames = true;

            % Import the data
            diagnos = readtable(diag_file, opts);

        catch
            warning(['file "',diag_file,'" incorrectly configured.'])
            continue
        end
        
        % find indices to keep
        if strcmp(files_to_clean{nn}, 'plot_times')
            time = diagnos.SimulationTime_s_;
        else
            time = diagnos.Time;
        end
        %%
        [~, uniqs] = unique(time);
        uniqlogical = zeros(size(time));
        uniqlogical(uniqs) = 1;
        findinds = zeros(size(time));
        findinds(find_inds(time)) = 1;
        %%
        keep_inds = findinds' & (isfinite(time))' & (uniqlogical)';
        %%
        
        % clean up diagnostics file
        names = fieldnames(diagnos);
        N_fields = length(names) - 1;
        new_diagnos = struct;
        for ii = 1:N_fields
            if ~isempty(strfind('Properties, Row', names{ii}))
                continue
            end
            temp = diagnos.(names{ii});
            new_diagnos.(names{ii}) = temp(keep_inds);
        end
        new_diagnos.Properties = diagnos.Properties;
        
        save([files_to_clean{nn},'.mat'],'-struct','new_diagnos')
        disp(['Completed cleaning of ',diag_file])
    else
        disp([diag_file,' not found']);
    end
end

end % of function

function keep_inds = find_inds(time)
% simplify time variable and find restarts
step_time = [time(1); time(2:end)-time(1:end-1)];
restart_inds = find(step_time<0);
last_inds = [restart_inds-1; length(time)];   % last indices before a restart
restart_time = time(restart_inds);
N_restarts = length(last_inds);

% find indices to keep (remove unwanted doubled time)
for ii = 1:N_restarts
    if ii == 1
        keep_inds = 1:last_inds(ii);
    elseif ii > 1
        old_inds = keep_inds(time(keep_inds) < restart_time(ii-1));
        add_inds = (last_inds(ii-1)+1):last_inds(ii);
        keep_inds = [old_inds add_inds];
    end
end
%keep_inds = keep_inds';
end % end of function
