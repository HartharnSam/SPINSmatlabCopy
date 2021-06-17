function rho = nleos(s,t,varargin)

%
% NLEOS.M    Density of sea water - adapted from densmdjwf.m for MITGCM
%=========================================================================
%
% USAGE:  dens = nleos(S,T,P)
%
% DESCRIPTION:
%    Density of Sea Water using the McDougall et al. 2003 (JAOT 20)
%    polynomial.
%
% INPUT:  (all must have same dimensions)
%   S     = salinity    [psu      (PSS-78)]
%   Theta = potential temperature [degree C (IPTS-68)]
%   P     = pressure    [dbar] (optional argument)
%   
%   Input variables should be the same size or scalars. 
%   In no P is supplied, assume sea surface pressure (p=0)
%
% OUTPUT:
%   dens = density  [kg/m^3] 
% 
% AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu)
%
% check value
% S     = 35.5 PSU
% Theta = 3 degC
% P     = 3000 dbar
% rho   = 1041.83267 kg/m^3


% McDougall et al., 2003, JAOT 20(5), pp. 730-741

% created by mlosch on 2002-08-09
% MODIFIED FOR UW FLUIDS GROUP: Andrew Grace 2021-05-31
%----------------------
% CHECK INPUT ARGUMENTS
%----------------------

  if nargin < 2
    error('nleos.m: Must pass at least 2 parameters, S and T.')
  elseif nargin > 3
    error('nleos.m: Must pass no more than 3 parameters, S and T and p.')
  end 

  %The following block of code accounts for when you have a temperature
  %(salinity) array, and scalar salinity (temperature). Scalar variable is changed
  %to size of array variable.
  
    [ms,ns] = size(s);
    [mt,nt] = size(t);

    if (ms==1 & ns==1) & (mt~=1 | nt~=1)  % S is a scalar and T is not.  
        s = s*ones(mt,nt); %Fill S to size of T
    elseif (mt==1 & nt==1) & (ms~=1 | ns~=1) % T is a scalar and T is not.
        t = t*ones(ms,ns); %Fill T to size of S
    elseif mt==ms & nt==ns     % T is a array of size(S)
       % shape ok 
    else
        error([' are incompatible sizes. Please enter the ' ...
               'same sized arrays, or one or both as scalar(s).'])
    end 

  if (length(varargin) == 1)
      p = varargin{1}; %If user enters pressure, assign it to p
      [mp,np] = size(p);
  else
      p = 0; %if user does not enter a pressure, assume sea surface. 
      [mp,np] = size(p);
  end
    if (mp==1 & np==1) & (mt~=1 | nt~=1)  % P is a scalar.  
        p = p*ones(mt,nt); %Fill P to size of T
    elseif mp==mt & np==nt     % P is a array of size(T) 
       %(Could have put S, but they are the same size at this point)
       % shape ok, do nothing.
    else
        error(['P is compatible with S and T. Please enter an array of size(S) ' ...
               ' os size(T), or as a scalar.'])
    end 

  % coefficients nonlinear equation of state in pressure coordinates for
  eosMDJWFnum(12) =  9.99843699e+02;
  eosMDJWFnum( 1) =  7.35212840e+00;
  eosMDJWFnum( 2) = -5.45928211e-02;
  eosMDJWFnum( 3) =  3.98476704e-04;
  eosMDJWFnum( 4) =  2.96938239e+00;
  eosMDJWFnum( 5) = -7.23268813e-03;
  eosMDJWFnum( 6) =  2.12382341e-03;
  eosMDJWFnum( 7) =  1.04004591e-02;
  eosMDJWFnum( 8) =  1.03970529e-07;
  eosMDJWFnum( 9) =  5.18761880e-06;
  eosMDJWFnum(10) = -3.24041825e-08;
  eosMDJWFnum(11) = -1.23869360e-11;
  
  eosMDJWFden(13) =  1.00000000e+00;
  eosMDJWFden( 1) =  7.28606739e-03;
  eosMDJWFden( 2) = -4.60835542e-05; 
  eosMDJWFden( 3) =  3.68390573e-07;
  eosMDJWFden( 4) =  1.80809186e-10;
  eosMDJWFden( 5) =  2.14691708e-03;
  eosMDJWFden( 6) = -9.27062484e-06;
  eosMDJWFden( 7) = -1.78343643e-10;
  eosMDJWFden( 8) =  4.76534122e-06;
  eosMDJWFden( 9) =  1.63410736e-09;
  eosMDJWFden(10) =  5.30848875e-06;
  eosMDJWFden(11) = -3.03175128e-16;
  eosMDJWFden(12) = -1.27934137e-17;

  p1 = p;
    
  t1 = t;
  t2 = t.*t;
  
  s1=s;
  is = find(s1(:) < 0 );
  if ~isempty(is)
    warning('found negative salinity values, reset them to NaN');
    s1(is) = NaN;
  end
  sp5 = sqrt(s1);
  p1t1=p1.*t1;
  %
  num = eosMDJWFnum(12) ...
	+ t1.*(eosMDJWFnum(1) ... 
        +     t1.*(eosMDJWFnum(2) + eosMDJWFnum(3).*t1) )  ...
	+ s1.*(eosMDJWFnum(4) ...
        +     eosMDJWFnum(5).*t1  + eosMDJWFnum(6).*s1) ...
	+ p1.*(eosMDJWFnum(7) + eosMDJWFnum(8).*t2 ...
        +     eosMDJWFnum(9).*s1 ...
	+     p1.*(eosMDJWFnum(10) + eosMDJWFnum(11).*t2) );
  den = eosMDJWFden(13) ...
	+ t1.*(eosMDJWFden(1) ...
	+     t1.*(eosMDJWFden(2) ...
	+         t1.*(eosMDJWFden(3) + t1.*eosMDJWFden(4) ) ) ) ...
	+ s1.*(eosMDJWFden(5) ...
	+     t1.*(eosMDJWFden(6) ... 
	+         eosMDJWFden(7).*t2) ... 
	+     sp5.*(eosMDJWFden(8) + eosMDJWFden(9).*t2) ) ... 
	+ p1.*(eosMDJWFden(10) ...
	+     p1t1.*(eosMDJWFden(11).*t2 + eosMDJWFden(12).*p1) );
  
  epsln = 0;
  denom = 1.0./(epsln+den) ;

  
  rho = num.*denom;

  return
  
