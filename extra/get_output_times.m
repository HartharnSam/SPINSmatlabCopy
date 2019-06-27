function [times, outputs] = get_output_times()
% Find the outputs and times associated with a particular run
%
% David Deepwell, 2019

id = 'MATLAB:table:ModifiedAndSavedVarnames';
warning('off',id)
pt = readtable('plot_times.txt');
warning('on',id)

times   = pt.SimulationTime_s_;
outputs = pt.OutputNumber;

% check for first output
if first_output() == outputs(1) - 1
    outputs = [first_output(); outputs];

    par = spins_params();
    t_0 = times(1) - par.plot_interval;
    times = [t_0; times];
elseif first_output() ~= outputs(1)
    warning('Something has gone wrong with time.')
end
