function rho = quadeos(t)

%
% QUADEOS.M    Density of cold freshwater
%=========================================================================
%
% USAGE:  dens = quadeos(T)
%
% DESCRIPTION:
%   Density of cold fresh water using a quadratic fit 
%   to the McDougall et al. 2003 (JAOT 20) polynomial equation of state.
%   Follows the form:
%   
%   rho(T) = rho(T_md) + C*(T - Tmd)^2
%   Note that this EOS should should not be used far away from the
%   Temperature of maximum density. When in doubt, check against nleos.m
%   with 0 salinity and 0 pressure.
%
% INPUT: 
%   T = in situ temperature [degree C]
%
% OUTPUT:
%   dens = density  [kg/m^3] 
% 
% AUTHOR:  Andrew Grace 2021-06-01  (andrew.grace@uwaterloo.ca)

if any(t>10)
    warning(['Not recommended for use for temperatures greater than 10 degrees C']);
end

C = -0.007641729395312;
rhomax = 9.999744074665388e+02;
Tmd = 3.973973973973974;

rho = rhomax + C*(t - Tmd).^2;
return 
