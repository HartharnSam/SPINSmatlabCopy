function change_dt(input_times)

fprintf('This script will delete data \n');
isContinue = input('Are you sure you wish to continue? Y/N [N]: ', 's');


if isempty(isContinue) || strcmpi(isContinue, 'n')
    fprintf('Script terminated \n')
    return
end

%% Now for the actual script
[~, outputs] = get_output_times(false);
keep_inds = clean_times(input_times);

remove_outputs = outputs(~keep_inds);
mkdir('archive')
for ii = 1:length(remove_outputs)
   movefile(['*.', num2str(remove_outputs(ii))], 'archive/')
end
end


function keep_inds = clean_times(input_times)
% clean_diagnostics reads in the diagnostic file and removed repeated
% times, and saves clean version as a .mat file

files_to_clean = 'plot_times'; % without file extensions (all are .txt)
%%
% read analysis file (if it exists
diag_file = [files_to_clean,'.txt'];
if exist(diag_file, 'file') == 2
    try
        diagnos = readtable(diag_file);
    catch
        warning(['file "',diag_file,'" incorrectly configured.'])
        
    end
    
    % find indices to keep
    time = diagnos.SimulationTime_s_;
    time = round(time, 1);
    %%
    keep_inds = any(input_times == time,2);

    % clean up diagnostics file
    new_diagnos = diagnos;
    new_diagnos(~keep_inds, :) = [];
    % Write to file
    writetable(new_diagnos, diag_file)
    disp(['Completed cleaning of ',diag_file])
else
    disp([diag_file,' not found']);
end


end % of function