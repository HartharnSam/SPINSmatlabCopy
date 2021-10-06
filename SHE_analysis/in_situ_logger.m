function [t, density, rho_normalised] = in_situ_logger(directory, x_loc)
%IN_SITU_LOGGER - simulates an in-situ density logger on the bed at point 
% x_loc (m)
%
% 
% Other m-files required: spinsgrid2d, nearest_index, spins_reader_new,
% get_output_times, spins_params
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 22-Sep-2021; Last revision: 22-Sep-2021
% MATLAB Version: 9.10.0.1739362 (R2021a) Update 5

cd(directory);
[x , ~] = spinsgrid2d;
x_ind = nearest_index(x(:, 1), x_loc);
[t, outputs] = get_output_times;
density = NaN(length(t), length(x_loc));
for i = 1:length(t)
    ti = outputs(i);
density(i, :) = spins_reader_new('rho', ti, x_ind, 1);
end
params = spins_params;
rho_normalised = (density+params.delta_rho/2)/params.delta_rho;
end