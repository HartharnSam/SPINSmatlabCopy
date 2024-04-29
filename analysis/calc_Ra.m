function Ra = calc_Ra()
%CALC_RA - Calculate the Rayleigh number for a SPINS Simulation
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    input1 - Description
%    input2 - Description
%    input3 - Description
%
% Outputs:
%    output1 - Description
%    output2 - Description
%
% Example:
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 12-Dec-2023; Last revision: 12-Dec-2023
% MATLAB Version: 9.10.0.1739362 (R2021a) Update 5

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------

params = spins_params;

if isfield(params, 'delta_rho')
    delta_rho = params.delta_rho;
elseif isfield(params, 'rho_top')
    delta_rho = abs(params.rho_top-params.rho_bot); 
else
    rho = spins_reader_new('rho', 0);
    delta_rho = max(rho(:))-min(rho(:));
end
L = params.h_pyc;
g = params.g;
mu = params.visco; 
kappa = params.kappa;

Ra = (delta_rho*L*L*L*g)./(mu*kappa);

%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------