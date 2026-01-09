function diff_plot(ii, xlimits, zlimits)
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
if nargin<2 || isempty(xlimits)
   xlimits = [0 params.Lx]+params.min_x;
end
if nargin<3
   zlimits = [0 params.Lz]+params.min_z;
end

if params.Ny == 1
    %% Set up x and y grids
    [x, z] = spinsgrid2d;
    %x = x-params.L_adj;

    %% Read in Data
    try
        rho = spins_reader_new('rho',ii);
    catch
        s = spins_reader_new('s', ii); t = spins_reader_new('t', ii);
        rho = nleos(s, t);
    end
    u = spins_reader_new('u',ii);
    umaxabs = max(abs(u(:)));

else
    [~, ~, ~, x, z] = spinsgrid3d;
    %x = x-params.L_adj;
    

    %% Read in Data
    rho = spins_reader_new('rho',ii, [], 1, []);
    u = spins_reader_new('u',ii, [], 1, []);
    umaxabs = max(abs(u(:)));
end
%% Set up figure
clf

%% Plot density first
s = subaxis(2,1,1);
pcolor(x,z,rho-mean(rho)),shading flat; %caxis(rhoRange);
hold on
rho_scale = (mean(rho)/max(abs(rho(:)))+1)*params.Lx/2;
plot(rho_scale, mean(z), 'w-');

title('Density (upper panel) and u (lower panel)')
%contourf(x,z,rho),shading flat % This is better for printed figures
ylabel('z (m)')
set(gca, 'xticklabels', {});
set(gca, 'xlim', xlimits);
set(gca, 'ylim', zlimits);

colormap(s, cmocean('dense'))
%colormap(s, 'jet')

%% Plot Horizontal Velocity Next
s = subaxis(2,1,2);
pcolor(x,z,u-mean(u)),shading flat,clim([-1 1]*umaxabs*0.05)
hold on
u_scale = (mean(u)/max(abs(u(:)))+1)*params.Lx/2;
plot(u_scale, mean(z), 'k-');
hold off
ylabel('z (m)')
xlabel('x (m)')
set(gca, 'xlim', xlimits);
set(gca, 'ylim', zlimits);

colormap(s, cmocean('balance'));


%% END OF CODE %%
% --------------------------------------------------