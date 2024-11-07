function [new_diagnos, movelist] = change_dt(input_times, isTest)
%CHANGE_DT - Reduces output frequency of SPINS outputs by updating
%plot_times.txt, moving unwanted outputs to archive, and renaming the
%existing outputs
% Has the same effect as if plot_interval had been set differently.
%
% Inputs:
%    input_times - times (not outputs) at which we want data kept for
%
% Other m-files required: get_output_times
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% Department of Geography & Environmental Sciences, Northumbria University
% email address: sam.hartharn-evans@northumbria.ac.uk
% GitHub: https://github.com/HartharnSam
% 28-Oct-2024; Last revision: 28-Oct-2024
% MATLAB Version: 23.2.0.2668659 (R2023b) Update 9

arguments
    input_times (1, :) double
    isTest (1, 1) logical
end

%% DOUBLE CHECK WE ACTUALLY WANT TO DO THE DAMAGE
if ~isTest
    fprintf('This script will delete data \n');
    isContinue = input('Are you sure you wish to continue? Y/N [N]: ', 's');

    if isempty(isContinue) || strcmpi(isContinue, 'n')
        fprintf('Script terminated, no data adjusted \n')
        return
    end
end

%% First, update plot_times.txt, and identify how to adjust outputs
% This part is based on clean_diagnostics.m - it reads in the plot_times.txt file and removes repeated
% times, and the ones we don't want and saves clean version.

if ~isTest
    mkdir('archive')
end

files_to_clean = 'plot_times'; % without file extensions (all are .txt)
clean_diagnostics();
% read analysis file (if it exists)
diag_file = [files_to_clean,'.txt'];
if exist(diag_file, 'file') == 2
    try
        diagnos = readtable(diag_file);
    catch
        error(['file "',diag_file,'" incorrectly configured.'])
    end

    % find indices to keep
    time = diagnos.SimulationTime_s_;
    og_outputs = diagnos.OutputNumber;
    % Check if there is a 0 time in the plot_times.txt file, add one if not
    add_zerotime = false;
    if (first_output() == og_outputs(1) - 1)
        og_outputs = [first_output(); og_outputs];
        par = spins_params();
        t_0 = max(0, time(1) - par.plot_interval);
        time = [t_0; time];
        add_zerotime = true;
    elseif (first_output() ~= og_outputs(1))
        warn_msg = sprintf(['The first output file does not match\n',...
            '    with the first output as stated in plot_times.txt']);
        warning(warn_msg);
    end

    % Round times for matching, and then check what indices are wanted in
    % the new regime.
    time = round(time, 1);
    input_times = round(input_times, 1);
    keep_inds = any(input_times == time,2);
    remapped_outputs = og_outputs(keep_inds);

    % clean up  file
    new_diagnos = diagnos;
    if add_zerotime
        new_diagnos(~keep_inds(2:end), :) = [];
    else
        new_diagnos(~keep_inds, :) = [];
    end
    new_diagnos.OutputNumber = (1:height(new_diagnos))';
    % Write to file
    if ~isTest
        writetable(new_diagnos, diag_file)
        copyfile(diag_file, "archive\");
    end
    disp(['Completed cleaning of ',diag_file])
else
    error([diag_file,' not found']);
end

%% First, remove the old outputs
remove_outputs = og_outputs(~keep_inds);
if ~isTest
    for ii = 1:length(remove_outputs)
        movefile(['*.', num2str(remove_outputs(ii))], 'archive/')
    end
end
disp("Moved old files");
%% Now re-map the old data
% First identify what old data needs remapping
output_0 = first_output();
list_of_outputs = dir("./*."+output_0);
list_of_outputs = extractfield(list_of_outputs, "name");
for jj = 1:length(list_of_outputs)
    tmp = strsplit(list_of_outputs{jj}, '.');
    output_list{jj} = tmp{1};
end
movelist = "";
for ii = 1:length(remapped_outputs)
    if ~isequal(remapped_outputs(ii), ii-1)
        for jj = 1:length(output_list)
            if isTest
                movelist(end+1, :) = "./"+output_list{jj}+"."+remapped_outputs(ii)+" to "+"./"+output_list{jj}+"."+(ii-1);
                if ~exist("./"+output_list{jj}+"."+remapped_outputs(ii), "file")
                    error("okay")
                end
            else
                try
                    movefile("./"+output_list{jj}+"."+remapped_outputs(ii), "./"+output_list{jj}+"."+(ii-1))
                    movelist(end+1, :) = "./"+output_list{jj}+"."+remapped_outputs(ii)+" to "+"./"+output_list{jj}+"."+(ii-1);
                catch
                    warning("Output "+(ii-1)+" check again - was "+remapped_outputs(ii));
                end
            end
        end
    end
end
%% As a sanity check, its worth deleting some common derivatives too, in case they arise at a later timestep than 0

delete vorty.* diss.* *_z.* *_x.* N2.* ri.*
