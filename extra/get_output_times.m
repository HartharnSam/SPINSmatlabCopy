function [times, outputs] = get_output_times(do_check)
% Find the outputs and times associated with a particular run
%
% David Deepwell, 2019

if nargin == 0
    do_check = true;
end

id = 'MATLAB:table:ModifiedAndSavedVarnames';
warning('off',id)
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["SimulationTime_s_", "OutputNumber", "VarName3", "VarName4"];
opts.VariableTypes = ["double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
pt = readtable('plot_times.txt', opts);
warning('on',id)

times   = pt.SimulationTime_s_;
inds = isnan(times);
times = times(~inds);
outputs = pt.OutputNumber;
outputs = outputs(~inds);
% check for first output
if do_check
    if first_output() == outputs(1) - 1
        outputs = [first_output(); outputs];

        par = spins_params();
        t_0 = times(1) - par.plot_interval;
        times = [t_0; times];
    elseif first_output() ~= outputs(1)
        warn_msg = sprintf(['The first output file does not match\n',...
            '    with the first output as stated in plot_times.txt']);
        warning(warn_msg);
    end
end
