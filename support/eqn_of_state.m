function rho = eqn_of_state(t,s)
% Computes density of seawater at 0 pressure by the equation
% of state in Brydon et al (1999)
%
% courtesy Fluids Group UW

if (~exist('s')) s=0; end

% Constants
c1 = -9.20601e-2;
c2 =  5.10768e-2;
c3 =  8.05999e-1;
c4 = -7.40849e-3;
c5 = -3.01036e-3;
c6 =  3.32267e-5;
c7 =  3.21931e-5;

rho = c1 + c2.*t + c3.*s + c4.*t.^2 + ...
      c5.*s.*t + c6.*t.^3 + c7.*s.*t.^2;
