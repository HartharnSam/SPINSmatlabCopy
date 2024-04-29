function plumeCharacteristics = characterise_plume(startframe,endframe, zlims)
%CHARACTERISE_PLUME - Carries out calculations of plume parameters
% Calculates timeseries of plume depth, velocity, density; Richardson
% number,
%
% Syntax:  plumeCharacteristics = characterise_plume(startframe,endframe)
%
% Inputs:
%    startframe - Description
%    endframe - Description
%
% Outputs:
%    plumeCharacteristics - Description
%
% Other m-files required: clencurt, contourf_clrbar, figure_print_format,
% find_position, first_output, get_grad2, last_output, spins_params,
% spins_reader_new, spinsgrid2d
%
% Subfunctions: none
% MAT-files required: none
%
% See also: characterize_wave
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 22-Nov-2023; Last revision: 22-Nov-2023
% MATLAB Version: 9.10.0.1739362 (R2021a) Update 5

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
close all;
params = spins_params;
%% Spatial grids
[x, z] = spinsgrid2d;
Nx = params.Nx; Nz = params.Nz;
dx = diff(x, 1, 1); dx(end+1, :) = dx(end);

if nargin <3 || isempty(zlims)
    zlims = [0 params.Lz]+params.min_z; % For hovmollers only
end

%% Time settings
% Find out which outputs exist
[time, og_outputs]= get_output_times(false);
first_out_full = og_outputs(1);
last_out_full = og_outputs(end);

if nargin == 0 || (isempty(startframe) && isempty(endframe))
    startframe = first_out_full;
    endframe = last_out_full;
end

% % Finalise timey things
first_out = max(first_out_full, startframe);
last_out = min(last_out_full, endframe);

output_inds = find(og_outputs == first_out):find(og_outputs == last_out);

time = time(output_inds);
outputs = og_outputs(output_inds);

noutputs = length(outputs);

%% Other settings
% Set saving filename
filename = 'plume_characteristics';

% Set/find isopycnal parameters
% Settings for one isopycnal
rho_1 = params.rho_bot; rho_2 = params.rho_top;

contval = rho_1 + 0.5*(rho_2 - rho_1); % Set as the mid-density currently, but open to change the 0.5 to any ratio

% Set up vectors of wave characteristics
n_cont = length(contval);
h_plume = NaN(last_out_full, n_cont);
top_rhobar = h_plume; bot_rhobar = h_plume;
top_ubar = h_plume; bot_ubar = h_plume;
Ri_prof = NaN(last_out_full, Nz);
param_times = NaN(last_out_full, 1);

rho_ts = NaN(size(z, 2), length(time));
ri_ts = rho_ts; u_ts = rho_ts; drhodz_ts = rho_ts; vorty_ts = rho_ts;

% Calculate Chebyschev volumes
% Compute the area associated with each Chebyshev point using the values
% halfway between the point below and above
Nzc = Nz-1;
%[~,z1dc] = cheb(Nzc);
[~,wci] = clencurt(Nzc);

% Then normalise this by the change in depth at each x position
winow = NaN(Nx, Nz);
for jj = 1:Nx
    Lznow = max(z(jj, :))-min(z(jj, :));
    %arcphys(jj, :) = arc*Lznow; % Thickness of each grid point
    % get the local chain rule expression - Weight of each grid point
    winow(jj, :) = wci*(Lznow)*0.5*dx(jj, 1);
end
% Vol = sum(winow(:)); % Sanity check for volume calculations


% Loop through outputs
%full_output_inds =

% pre-slicing for parfor
Ny = params.Ny; Lz = params.Lz; min_z = params.min_z; g = params.g;
rho_0 = params.rho_0; slope = params.slope;
z_mid = z(Nx/2,:);

% Start parallel pool things
pool = gcp;
D = parallel.pool.DataQueue;
h = waitbar(0, 'Characterising ...');
afterEach(D, @nUpdateWaitbar);

p = 1; N = noutputs;
    function nUpdateWaitbar(~)
        waitbar(p/N, h);
        p = p + 1;
    end

warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');

parfor jj = 1:noutputs
    ii = outputs(jj);
    current_time = time(jj);
    full_output_ind(jj) = find(og_outputs == ii);
    % Read data and find reference stratification
    if Ny > 1 % in 3 dimensions
        warning('not yet configured properly')
        %rho = spins_reader_new('Mean Density', ii, x_inds, 1, z_inds);
        %u = spins_reader_new('Mean u', ii, x_inds, 1, z_inds);
    else
        rho = spins_reader_new('rho', ii);
        u = spins_reader_new('u', ii);
    end
    strat = mean(rho);
    
    for nn = 1:n_cont
        % find background depth of chosen isopycnal (contval)
        [strat_pos, ~] = find_position(z_mid, strat, contval(nn));         %#ok<PFBNS>
        h_plume_tmp(jj, nn) = Lz+min_z - strat_pos;
        %[cont_x, cont_y] = find_contour(x, z, rho, contval(nn));
        
        top_inds = rho < contval(nn);
        bot_inds = rho > contval(nn);
        
        top_rhobar_tmp(jj, nn) = sum((rho(top_inds).*winow(top_inds)))/sum(winow(top_inds), 'all'); %#ok<PFBNS>
        bot_rhobar_tmp(jj, nn) = sum((rho(bot_inds).*winow(bot_inds)))/sum(winow(bot_inds), 'all');
        
        top_ubar_tmp(jj, nn) = sum((u(top_inds).*winow(top_inds)))/sum(winow(top_inds), 'all');
        bot_ubar_tmp(jj, nn) = sum((u(bot_inds).*winow(bot_inds)))/sum(winow(bot_inds), 'all');
        
    end
    % Add in a mean layer density measure
    % And a layer velocity measure
    
    % set-up the figure
    %all_conts = true; % do all contours exist in the domain?
    clf
    hold on
    
    g_rho0 = -g/rho_0;
    [~, du_dz] = get_grad2(u);
    [~, drho_dz] = get_grad2(rho);
    N_sq = g_rho0 * drho_dz.*cosd(slope);
    Ri = N_sq./(du_dz.^2);
    Ri_prof_tmp(jj, :) = mean(Ri);
    % TODO: Add the things L178:278 in characterize_wave
    param_times(jj) = current_time;
    %completion(jj, noutputs, .1, 'Characterising Plume');
    
    % Calculate mean profiles
    rho_ts(:, jj) = strat;
    u_ts(:, jj) = mean(u);
    vorty_ts(:, jj) = mean(spins_derivs('vorty', ii));
    ri_ts(:, jj) = Ri_prof_tmp(jj, :);
    drhodz_ts(:, jj) = mean(drho_dz);
    
    send(D, jj);
end

h_plume(full_output_ind, :) = h_plume_tmp(1:noutputs, :);
top_rhobar(full_output_ind, :) = top_rhobar_tmp(1:noutputs, :);
bot_rhobar(full_output_ind, :) = bot_rhobar_tmp(1:noutputs, :);
top_ubar(full_output_ind, :) = top_ubar_tmp(1:noutputs, :);
bot_ubar(full_output_ind, :) = bot_ubar_tmp(1:noutputs, :);
Ri_prof(full_output_ind, :) = Ri_prof_tmp(1:noutputs, :);

Ri_b = -9.81*(cosd(params.slope)*h_plume.*(top_rhobar-bot_rhobar))./((top_ubar-bot_ubar).^2);

%% Plot the characteristics plots
first_time = time(1); last_time = time(end);

figure('Name', 'Timeseries');
tiledlayout('flow', 'TileSpacing','tight', 'Padding','tight');
nexttile;
plot(param_times, h_plume);
xlabel('time (s)'); ylabel('h (m)');
xlim([first_time last_time])

nexttile;
plot(param_times, (top_rhobar-top_rhobar(1))*1000); hold on
plot(param_times, -bot_rhobar*1000);
xlabel('time (s)'); ylabel('$\overline{{\rho}}$');
axis padded
xlim([first_time last_time])

nexttile;
plot(param_times, top_ubar); hold on
plot(param_times, -bot_ubar);
xlabel('time (s)'); ylabel('$\overline{{u}}$');
xlim([first_time last_time])

nexttile;
plot(param_times, Ri_b);
hold on
yline(0.25);
xlabel('time (s)'); ylabel('$Ri$');
ylim([0 1]);
xlim([first_time last_time])

nexttile; nexttile;
inds = ~isnan(param_times);
cc = contourf_clrbar(param_times(inds), z(1, :), Ri_prof(inds, :)', [-Inf 0 .25 1 2 Inf], 'cmocean(amp');
cc.LineColor = 'none';
figure_print_format(gcf, 12);
exportgraphics(gcf, [filename, '.png']);

%% Calculate hovmoller plots
rizlims = max(z(:))- h_plume*[3 0];
zprof = z(1, :);

figure('Name', 'Hovmollers');
tiledlayout(6, 1);

% Plot mean density profile
nexttile;
pcolor(param_times, zprof, rho_ts);
hold on
plot(param_times, params.Lz-h_plume, ':k');
ylim(zlims);
c = colorbar;
ylabel(c, '$\rho$');
ylabel('$z (m)$');
caxis([params.rho_top params.rho_bot]);
cmocean('dense');

% Plot density gradient
nexttile
pcolor(param_times, zprof, log10(drhodz_ts.^2)); %TODO: Should this be a log?
ylim(zlims);
c = colorbar;
caxis([-10 0]);
ylabel(c, '$\frac{d\rho}{dz}$');
ylabel('$z (m)$')
cmocean('thermal');

nexttile;
pcolor(param_times, zprof, u_ts);
ylim(zlims);
c = colorbar;
caxis([-.5 .5]);
ylabel(c, '$u (m/s)$');
ylabel('$z (m)$')
cmocean('balance', 'pivot', 0);

nexttile;
pcolor(param_times, zprof, vorty_ts);
ylim(zlims);
c = colorbar;
ylabel(c, '$\omega (s^{-1})$');
ylabel('$z (m)$')
cmocean('balance', 'pivot', 0);

nexttile;
pcolor(param_times, zprof, ri_ts);
ylim(zlims); caxis([0 1])
ylabel('$z (m)$')
c = colorbar;
ylabel(c, '$Ri$');
cmocean('amp');
hold on
%yline(rizlims(1));
plot(param_times, rizlims(:, 1));
hold off

nexttile;
zind = nearest_index(zprof', rizlims(:, 1)');
ri_ts_max = param_times.*NaN;
for ii = 1:length(zind)
    ri_ts_max(ii) = max(ri_ts(zind(ii):end, ii));
end

plot(param_times, ri_ts_max);
ylim([0 1]);
hold on; yline(.25, ':');
xlim([param_times(1) param_times(end)]);
xlabel('t (s)')
ylabel('$Ri$')

figure_print_format(gcf, 12);
set(gcf, 'Position', [681 91 560 842]);
exportgraphics(gcf, ['hovmollers', '.png'], 'resolution', 300);
%% Calculate some other dimensionless numbers
%% Calculate and print some dimensionless numbers
g = params.g;
H = params.Lz;
nu = params.visco;
rho_i = params.rho_top;
rho_a = params.rho_bot;

Gr = abs(g*((rho_i - rho_a))*(H^3)/(nu^2));
u_1 = top_ubar(full_output_ind(end));
Re = abs(u_1).*H/nu;
Ar = Gr/(Re.^2);

fprintf('\n \n --------- \n');
fprintf('Grashof =  %.2e | ', Gr);
if Gr > 1e10
    fprintf('Probably turbulent \n');
elseif Gr < 1e8
    fprintf('Probably laminar \n');
else
    fprintf('Marginal turbulent/laminar \n');
end

fprintf('Reynolds =  %.2e \n', Re);
fprintf('Archimedes = %2.2f | ', Ar);
if Ar > 1.5
    fprintf('Natural Convection \n');
elseif Ar < 0.8
    fprintf('Forced Convection \n');
else
    fprintf('Marginally Forced/Natural convection \n');
end
%fprintf('Richardson =  %2.2f \n', Ri);

%% Now save
%tmp_time = NaN(noutputs+outputs(1)-1, 1);
%tmp_time(outputs) = time;
%time = tmp_time;

if nargout == 1
    plumeCharacteristics = load(filename);
end
time = param_times;
h_plume_prime = time.*NaN;
h_plume_prime(~isnan(time)) = FiniteDiff(time(~isnan(time)),1,2,false,false)*...
    h_plume((~isnan(time)));

save(filename, 'time', 'h_plume', 'h_plume_prime', 'top_rhobar', 'bot_rhobar', 'top_ubar', 'bot_ubar', 'Ri_b', 'Ar', 'Gr', 'Re');

end
%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------

