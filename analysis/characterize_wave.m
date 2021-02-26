function WaveStats = characterize_wave(isTwoLayer, time_window)

% CHARACTERIZE_WAVE   track the wave core assuming
%                     it moves from left to right
%
% Uses an isopycnal to find:
%    - the location of the wave
%        - from location of maximum vertical displacement
%    - the wave amplitude
%        - from the vertical isopycnal displacement
%    - the wavelength
%        - from the horizontal distance for the ispycnal to
%          drop to half-max value
%
% Inputs:
%    isTwoLayer - [OPTIONAL Boolean] - is wave a two layer case or not (affects choice
%    of isopycnal). Default false
%    time_window - [OPTIONAL 1x2 vector] - Frame numbers to characterize wave based on (e.g. [5
%    50]). Default [0 last_output]
%
% Outputs:
%    WaveStats - Structure containing key wave measurements (amplitude,
%    wavelength, speed, time it hits the slope)
%
% David Deepwell, 2017

%% Set-up
% find which outputs exist
first_out = first_output();
last_out = last_output();
init_outputs = first_out:last_out;
%%
if nargin <2
    startframe = 0;
    endframe = length(init_outputs)-1;
else
    startframe = time_window(1);
    endframe = time_window(2);
end
if nargin < 1
    isTwoLayer = false;
end

%%
first_out = max(first_output(), startframe);
last_out = min(last_out(), endframe);
outputs = first_out:last_out;
noutputs = length(outputs);
% set filename
filename = 'wave_characteristics';

% read grid and parameters
gd.z = zgrid_reader;
gd.x = xgrid_reader;
params = spins_params;
%gdpar = spins_gridparams('Full');
%split_gdpar
% shorten parameters
%noutputs = length(dir('u.*'));
gdnames = fieldnames(gd);
params.ndims = length(gdnames);

Nx = params.Nx;
Nz = params.Nz;

% find background stratification
if params.ndims == 3
    strat = mean(spins_reader('rho', 0, Nx, [], []));
elseif params.ndims == 2
    strat = spins_reader_new('rho', 0, Nx/2, 1, []);
end

% isopycnal parameters
%   One isopycnal:
contval = 0; %#ok
config = {'depr'};
if isTwoLayer
    isopyc_loc = params.pyc_loc -params.h_halfwidth + 0.01;
else
    isopyc_loc = params.pyc_loc -params.h_halfwidth; %+ 0.01;
    
end
if isvector(gd.z)
    contval = interp1(gd.z, strat, isopyc_loc);
else
    contval = interp1(gd.z(Nx/2,:), strat, isopyc_loc);
end

% contval = [-1 1]*val;
% config = {'depr'}; %{'elev', 'depr'}; % whether the isopycnals correspond
% to waves of elevation ('elev') or depression ('depr')
% Depression waves will be flipped to appear as elevation
% waves. This makes characterization easier.

% initial reading indices
xlind = 1;
xrind = Nx;
zbind = 1;
ztind = Nz;
mult = 2; % multiplier of wavelength to compute region to load

% set-up vectors of wave characteristics
n_cont = length(contval);
amplitude = zeros(noutputs, n_cont);
wave_center =      amplitude;
wavelength_right = amplitude;
wavelength_left =  amplitude;
strat_loc =        amplitude;
reach_end = false;

% colormap of plot
if n_cont <= 2
    c_map = [rgb('Black'); rgb('Red')];
else
    c_map = darkjet(n_cont);
end


% and the associated times
%time = get_output_times();
time = (startframe:endframe)';

%% Loop through outputs

for jj = 1:noutputs
    if reach_end
        disp(['Wave has reached the tank end, skipping output '...
            ,num2str(ii+1),' and all after'])
        completion(jj, noutputs)
        continue
    end
    %disp(num2str(jj));
    % current output and indices to use
    ii = outputs(jj);
    x_inds = xlind:xrind;
    z_inds = zbind:ztind;
    
    
    
    for nn = 1:n_cont
        % find background depth of chosen isopycnal (contval)
        if isvector(gd.z)
            [strat_pos, ~] = find_position(gd.z, strat, contval(nn));
        else
            [strat_pos, ~] = find_position(gd.z(Nx/2,:), strat, contval(nn));
        end
        strat_loc(jj, nn) = strat_pos;
    end
    
    % read data
    if params.ndims == 3
        rho = spins_reader_new('Mean Density', ii, x_inds, 1, z_inds);
    elseif params.ndims == 2
        rho = spins_reader_new('rho', ii, x_inds, 1, z_inds);
    end
    
    % set-up the figure
    all_conts = true; % do all contours exist in the domain?
    figure(1), clf
    hold on
    
    % loop through the contours
    p_hand = gobjects(1, length(n_cont));
    for nn = 1:n_cont
        if contval(nn) < max(rho(:)) && contval(nn) > min(rho(:))
            
            % find contour (isopycnal)
            if isvector(gd.x)
                [cont_x, cont_y] = find_contour(gd.x(x_inds), gd.z(z_inds), rho', contval(nn));
            else
                [cont_x, cont_y] = find_contour(gd.x(x_inds,z_inds), gd.z(x_inds,z_inds), rho, contval(nn));
            end
            if isempty(cont_x) 
                if jj > noutputs*.75
                    reach_end = true;
                end
            end
            cont_y = cont_y - strat_loc(jj, nn); % shift so that cont_y=0 at far field
            % if wave is a depression, then flip vertically
            if strcmp(config{nn}, 'depr')
                cont_y = -cont_y;
            end
            
            % find amplitudes and locations of local maxima
            [max_val, max_pos, max_ind] = find_wave_max(cont_x, cont_y);
            amplitude(jj, nn)   = max_val(1);
            wave_center(jj, nn) = max_pos(1);
            % what to do about other components?
            
            % find wavelengths
            inds = max_ind(1):length(cont_x);
            if length(inds) > 1 % if wave hasn't reach tank end
                % fore (right) wavelength
                [front_loc, ~] = find_position(cont_x(inds), cont_y(inds), max_val(1)/2);
                wavelength_right(jj, nn) = front_loc - max_pos(1);
                % aft (left) wavelength
                inds = 1:max_ind;
                [back_loc, ~] = find_half_max(cont_x(inds), cont_y(inds));
                wavelength_left(jj, nn) = max_pos(1) - back_loc;
            else % wave has reached tank end
                wavelength_right(jj, nn) = NaN;
                wavelength_left(jj, nn)  = NaN;
                if jj > noutputs*.75
                reach_end = true;
                end
            end
            
            if isempty(cont_x)
                continue
            end
            % make plot to check
            p_hand(nn) = plot(cont_x, cont_y, '.', 'Color', c_map(nn,:));
            % plot wave centre and wavelengths
            if ~isnan(max_pos)
                plot([1 1]*max_pos(1), [0 max_val(1)],'-', 'Color', c_map(nn,:))
                plot([1 1]*front_loc,  [0 max_val(1)],'-', 'Color', c_map(nn,:))
                plot([1 1]*back_loc,   [0 max_val(1)],'-', 'Color', c_map(nn,:))
            end
        else
            warning(['contour ',num2str(contval(nn)),' is not in the data set, skipping.'])
            all_conts = false;
        end
    end
    % other stuff
    if isempty(cont_x)
        continue
    end
    grid on
    title(['t=',num2str(time(jj)),' s'])
    if all_conts 
        leg = arrayfun(@(contval) ['rho = ',num2str(contval)],contval,'Uni',0);
        legend(p_hand, leg)
    end
    drawnow
    hold off
    
    % update reading indices
    if jj > 1 && ~isnan(wave_center(jj, 1))
        travel = max((wave_center(jj, :) - wave_center(jj-1, :)));
        xr = max(wave_center(jj, :)) + travel + mult*max(wavelength_right(jj,:));
        xl = min(wave_center(jj, :))      - 0.5*mult*max(wavelength_left(jj, :));
        %zt = 2*max_val(1) + strat_pos;
        if isvector(gd.x)
            xlind = nearest_index(gd.x, xl);
            xrind = nearest_index(gd.x, xr);
            %ztind = nearest_index(gd.z, zt);
        else
            xlind = nearest_index(gd.x(:,1), xl);
            xrind = nearest_index(gd.x(:,1), xr);
            %ztind = nearest_index(gd.z(Nz,:), zt);
        end
        %zbind = 1;
    else
        xlind = 1;
        xrind = Nx;
    end
    
    % print percentage of completion
    completion(jj, noutputs)
end

% get other important information
wave_speed = 0*amplitude;
for nn = 1:n_cont
    wave_speed(:,nn) = FiniteDiff(time,1,2,false,false)*wave_center(:,nn);
end

% calculate time means for key parameters - over flat topography
if isfield(params, 'hill_slope')
    end_of_slope = params.Lx-((params.hill_height/params.hill_slope)+params.hill_end_dist);
elseif isfield(params, 'ice_length')
    end_of_slope = params.Lx-params.ice_length;
end

read_inds = find(wave_center>(params.L_adj*1.15) & wave_center<end_of_slope);

mean_amp = mean(amplitude(read_inds));
mean_speed = mean(wave_speed(read_inds));
end_of_slope = time(max(read_inds));
mean_wavelength  = mean(wavelength_right(read_inds)) + mean(wavelength_left(read_inds));

% Display and save data
WaveStats = struct('meanAmp', mean_amp, 'meanWaveSpeed', mean_speed,...
    'endSlope', end_of_slope, 'meanWavelength', mean_wavelength);

disp(['Avg. Amplitude = ', num2str(mean_amp), 'm'])
disp(['Avg. Wave Speed = ', num2str(mean_speed), 'm/s']);
disp(['Avg. Wave Length = ', num2str(mean_wavelength), 'm']);
if max(read_inds)<length(wave_center)
    disp(['Wave reached slope at ', num2str(end_of_slope), 's']);
end

% save data
save(filename,'time','amplitude','wave_center',...
    'wavelength_right','wavelength_left','strat_loc','wave_speed','contval', 'isopyc_loc', 'WaveStats');






