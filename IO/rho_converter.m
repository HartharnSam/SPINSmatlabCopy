function rho_new = rho_converter(rho)
%RHO_CONVERTER - Converts from normalised density to real world density
%
% Syntax:  [rho_new] = rho_converter(rho)
%
% Inputs:
%    rho - Density field from SPINS (using spins_reader_new)
%
% Outputs:
%    rho_new - Adjusted density field from SPINS in kg/m^3
%
% Other m-files required: spins_params
% Subfunctions: none
% MAT-files required: none
%
% See also: SPINS_READER_NEW
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% June-2020; Last revision: 26-Nov-2020
% MATLAB Version: 9.9.0.1467703 (R2020b)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------

params = spins_params;
try
    delta_rho = params.delta_rho;
catch
    delta_rho = params.delta_rho_1+params.delta_rho_2;
end
rho_0 = params.rho_0;

rho_new = rho_0.*(1 + delta_rho/2 + rho.*(delta_rho.*(1+delta_rho) + 1));
end

