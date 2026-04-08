function times = outputs2times(outputs)
% Converts simulation outputs to times
% See also: times2outputs
% Sam Hartharn-Evans
[sim_times, sim_outputs] = get_output_times(true);

times = outputs.*NaN;
for ii = 1:length(outputs)
    times(ii) = sim_times(sim_outputs==outputs(ii));
end
