function [alpha,beta,rho_0] = compute_eos_coeffs(s0,t0)

%
% COMPUTE_EOS_COEFFS.M   Coefficients for linear equation of state
%=========================================================================
%
% USAGE:  [alpha,beta,rho_0] = compute_eos_coeffs(S0,T0)
%
% DESCRIPTION:           Computes reference density, thermal expansion
%                        coefficient and haline contraction coefficient 
%                        based on user supplied reference salinity and 
%                        temprature.
%
% INPUT: 
%   T0 = scalar reference temperature [degree C]
%   S0 = scalar reference Salinity [psu]
% OUTPUT:
%   alpha = thermal expansion coefficient  [kg/DegC.m^3] 
%   beta = haline contraction coefficient  [kg/m^3] 
%   rho_0 = reference density  [kg/m^3] 
%
% AUTHOR:  Andrew Grace 2021-06-01  (andrew.grace@uwaterloo.ca)



dT = 1e-10;
dS = 1e-10;

  
[ms,ns] = size(s0);
[mt,nt] = size(t0);

  if ms ~= 1 | ns ~=1 | mt ~=1 | nt ~=1
    error('compute_eos_coeffs.m: S0 and T0 must be scalars')
  end 


rho_0 = nleos(s0,t0);
alpha = (nleos(s0,t0 + dT) - nleos(s0,t0))/dT;
beta = (nleos(s0 + dS,t0) - nleos(s0,t0))/dS;

return 
