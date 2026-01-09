function outputs = times2outputs(times)
% Converts simulation times to output number
% Sam Hartharn-Evans
[sim_times, sim_outputs] = get_output_times(true);

outputs = times.*NaN;
for ii = 1:length(times)
    outputs(ii) = sim_outputs(nearest_index(sim_times, times(ii)));
end
