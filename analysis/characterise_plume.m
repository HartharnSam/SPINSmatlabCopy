function plumeCharacteristics = characterise_plume(startframe,endframe)
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
[x, z] = spinsgrid2d;
Nx = params.Nx; Nz = params.Nz;

dx = diff(x, 1, 1); dx(end+1, :) = dx(end);

%% Time settings
% Find out which outputs exist
[time, outputs]= get_output_times(false);
first_out_full = outputs(1);
last_out_full = outputs(end);

if nargin == 0
    startframe = first_out_full;
    endframe = last_out_full;
end

% % Finalise timey things
first_out = max(first_out_full, startframe);
last_out = min(last_out_full, endframe);

time = time((first_out:last_out)-outputs(1)+1);
outputs = first_out:last_out;

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
for jj = 1:noutputs
    ii = outputs(jj);
    current_time = time(jj);
    % Read data and find reference stratification
    if params.Ny > 1 % in 3 dimensions
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
        [strat_pos, ~] = find_position(z(Nx/2,:), strat, contval(nn));
        h_plume(ii, nn) = params.Lz+params.min_z - strat_pos;
        %[cont_x, cont_y] = find_contour(x, z, rho, contval(nn));
        top_inds = rho < contval(nn);
        bot_inds = rho > contval(nn);
        
        top_rhobar(ii, nn) = sum((rho(top_inds).*winow(top_inds)))/sum(winow(top_inds), 'all');
        bot_rhobar(ii, nn) = sum((rho(bot_inds).*winow(bot_inds)))/sum(winow(bot_inds), 'all');
        
        top_ubar(ii, nn) = sum((u(top_inds).*winow(top_inds)))/sum(winow(top_inds), 'all');
        bot_ubar(ii, nn) = sum((u(bot_inds).*winow(bot_inds)))/sum(winow(bot_inds), 'all');
        
    end
    % Add in a mean layer density measure
    % And a layer velocity measure
    
    % set-up the figure
    %all_conts = true; % do all contours exist in the domain?
    clf
    hold on
    
    g_rho0 = -params.g/params.rho_0;
    [~, du_dz] = get_grad2(u);
    [~, drho_dz] = get_grad2(rho);
    N_sq = g_rho0 * drho_dz.*cosd(params.slope);
    Ri = N_sq./(du_dz.^2);
    Ri_prof(ii, :) = mean(Ri);
    % TODO: Add the things L178:278 in characterize_wave
    param_times(ii) = current_time;
    completion(jj, noutputs, .1, 'Characterising Plume');
end

Ri_b = -9.81*(cosd(params.slope)*h_plume.*(top_rhobar-bot_rhobar))./((top_ubar-bot_ubar).^2);

%% Calculate some other dimensionless numbers
%% Calculate and print some dimensionless numbers
g = params.g;
H = params.Lz;
nu = params.visco;
rho_i = params.rho_top;
rho_a = params.rho_bot;

Gr = abs(g*((rho_i - rho_a))*(H^3)/(nu^2));
u_1 = top_ubar(end);
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


%% Now save and plot
%tmp_time = NaN(noutputs+outputs(1)-1, 1);
%tmp_time(outputs) = time; 
%time = tmp_time;

first_time = time(1); last_time = time(end);

figure('Name', 'Timeseries');
tiledlayout('flow');
nexttile;
plot(param_times, h_plume);
xlabel('time (s)'); ylabel('h (m)');
xlim([first_time last_time])

nexttile;
plot(param_times, top_rhobar*1000); hold on
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
xlabel('time (s)'); ylabel('$Ri$');
ylim([0 1]);
xlim([first_time last_time])

nexttile;
inds = ~isnan(param_times);
cc = contourf_clrbar(param_times(inds), z(1, :), Ri_prof(inds, :)', [-Inf 0 .25 1 2 Inf], 'cmocean(amp');
cc.LineColor = 'none';
figure_print_format(gcf, 14);
exportgraphics(gcf, [filename, '.png']);

if nargout == 1
    plumeCharacteristics = load(filename);
end
time = param_times;
h_plume_prime = time.*NaN; 
h_plume_prime(~isnan(time)) = FiniteDiff(time(~isnan(time)),1,2,false,false)*...
        h_plume((~isnan(time)));

save(filename, 'time', 'h_plume', 'h_plume_prime', 'top_rhobar', 'bot_rhobar', 'top_ubar', 'bot_ubar', 'Ri_b', 'Ar', 'Gr', 'Re');


%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------
