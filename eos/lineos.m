function rho = lineos(s,t,s0,t0)
%
% LINLEOS.M    Density of water (fresh or saline)
%=========================================================================
%
% USAGE:  rho = lineos(S,T,S0,T0)
%
% DESCRIPTION: Calculates rho as a linear function of temperature and
%              salinity of the form
%
%              rho = rho_0*(1 + alpha/rho_0*(T-T0) + beta/rho_0*(S-S0));
%
%              alpha, beta, and rho_0 are calculated by
%              compute_eos_coeffs.m
%
% INPUT: 
%   T = in situ temperature [degree C]
%   S = salinity
%   T0 = reference temperature
%   S0 = reference salinity
%
% OUTPUT:
%   dens = density  [kg/m^3] 
% 
% AUTHOR:  Andrew Grace 2021-06-01  (andrew.grace@uwaterloo.ca)

[ms,ns] = size(s0);
[mt,nt] = size(t0);

  if ms ~= 1 | ns ~=1 | mt ~=1 | nt ~=1
      error('lineos.m: S0 and T0 must be scalars')
  end 
  
[ms,ns] = size(s);
[mt,nt] = size(t);

    if (ms==1 & ns==1) & (mt~=1 | nt~=1)  % S is a scalar and T is not.  
        s = s*ones(mt,nt); %Fill S to size of T
    elseif (mt==1 & nt==1) & (ms~=1 | ns~=1) % T is a scalar and T is not.
        t = t*ones(ms,ns); %Fill T to size of S
    elseif mt==ms & nt==ns     % T is a array of size(S)
       % shape ok 
    else
        error(['lineos.m: S and T are incompatible sizes. Please enter the ' ...
               'same sized arrays, or one or both as scalar(s).'])
    end 

[alpha,beta,rho_0] = compute_eos_coeffs(s0,t0);

rho = rho_0*(1 + alpha/rho_0*(t - t0) + beta/rho_0*(s - s0));

return