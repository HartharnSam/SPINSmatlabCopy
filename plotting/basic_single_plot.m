function basic_single_plot(ii, xlimits)
%BASIC PLOTTING - produces simplest single frame plot of real density and horizontal
%velocity. 
%
% Inputs:
%    ii - frame number to be plotted [REQUIRED]
%    xlimits - X limits in real coordinates (e.g. [5 7]) [OPTIONAL]
%
% Other m-files required: subaxis
% Subfunctions: none
% MAT-files required: none
%
% See also: single_time_pic_example,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 05-Nov-2020; Last revision: 05-Nov-2020
% MATLAB Version: 9.9.0.1467703 (R2020b)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------

%% Read in parameters
params = spins_params;

%% Set up x and y grids
[x, z] = spinsgrid2d;
x = x-params.L_adj;
if nargin<2
    xlimits = [min(x(:, 1)) max(x(:, 1))];
end
    
%% Read in Data
rho=spins_reader_new('rho',ii);
rho = rho_converter(rho); % convert density to real density
rhoRange = [params.rho_0 params.rho_0+(params.delta_rho*params.rho_0)];
u = spins_reader_new('u',ii);
umaxabs = max(abs(u(:)));

%% Set up figure
figure
clf
colormap darkjet

%% Plot density first 
subaxis(2,1,1)
pcolor(x,z,rho),shading flat; caxis(rhoRange);
title('Density (upper panel) and u (lower panel)')
%contourf(x,z,rho),shading flat % This is better for printed figures
ylabel('z (m)')
set(gca, 'xticklabels', {});
set(gca, 'xlim', xlimits);

%% Plot Horizontal Velocity Next
subaxis(2,1,2)
pcolor(x,z,u),shading flat,caxis([-1 1]*umaxabs)
ylabel('z (m)')
xlabel('x (m)')
set(gca, 'xlim', xlimits);
%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------