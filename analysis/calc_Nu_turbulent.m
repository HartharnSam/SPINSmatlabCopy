function Nu = calc_Nu_turbulent(Ra, Pr)
%CALC_NU_TURBULENT - Calculates the Nusselt (or Sherwood) number for
%turbulent convection on a horizontal (or inclined) plate under turbulent
%conditions. 
% From Vliet, G.C., 1969, “Natural Convection Local Heat Transfer on 
% Constant-Heat-Flux Inclined Surface,” ASME J. Heat Transfer, Vol.
% 91, pp. 511-516.
%
% Inputs:
%    Ra - Rayleigh number (from compute_Ra)
%    Pr - Prandl Number (or Schmidt number to compute Sherwood number). If
%    blank, calculates based on available spins.conf files
%
% Outputs:
%    Nu - Nusselt Number (or Sherwood number if density)
%
% Other m-files required: spins_params
% Subfunctions: none
% MAT-files required: none
%
% See also: calc_Ra
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 12-Dec-2023; Last revision: 12-Dec-2023
% MATLAB Version: 9.10.0.1739362 (R2021a) Update 5


if nargin < 2
    params = spins_params;
    Pr = params.visco/params.kappa;
end

factor_1 = (0.387*(Ra^(1/6)));
factor_2 = (1+(0.492/Pr)^(9/16))^(8/27);

Nu = (0.825 + factor_1/factor_2)^2;

