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
% David Deepwell, 2017

%% Set-up

% set filename
filename = 'wave_characteristics';

% read grid and parameters
gdpar = spins_gridparams('Full');
split_gdpar
% shorten parameters
noutputs = params.noutputs;
Nx = params.Nx;
Nz = params.Nz;

% isopycnal parameters
%   One isopycnal:
%contval = 0;
%config = {'depr'};
%   Two isopycnals: (with finding density at a particular depth)
middepth = params.min_z + params.Lz/2;
h = params.h_halfwidth;
isopyc_loc = middepth - h;
strat = spins_reader('rho', 0, Nx, 1, []);
if isvector(gd.z)
    val = interp1(gd.z, strat, isopyc_loc);
else
    val = interp1(gd.z(Nx,:), strat, isopyc_loc);
end
contval = [-1 1]*val;
config = {'elev', 'depr'}; % whether the isopycnals correspond
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

% find which outputs exist
first_out = first_output();
last_out = last_output();
outputs = first_out:last_out;

%% Loop through outputs
for jj = 1:noutputs
    if reach_end
        disp(['Wave has reached the tank end, skipping output '...
        ,num2str(ii+1),' and all after'])
        completion(jj, noutputs)
        continue
    end

    % current output and indices to use
    ii = outputs(jj);
    x_inds = xlind:xrind;
    z_inds = zbind:ztind;

    % find background stratification
    if params.ndims == 3
        strat = mean(spins_reader('rho', ii, Nx, [], []));
    elseif params.ndims == 2
        strat = spins_reader('rho', ii, Nx, 1, []);
    end

    for nn = 1:n_cont
        % find background depth of chosen isopycnal (contval)
        if isvector(gd.z)
            [strat_pos, ~] = find_position(gd.z, strat, contval(nn));
        else
            [strat_pos, ~] = find_position(gd.z(Nx,:), strat, contval(nn));
        end
        strat_loc(jj, nn) = strat_pos;
    end

    % read data
    if params.ndims == 3
        rho = spins_readdata('Mean Density', ii, x_inds, 1, z_inds);
    elseif params.ndims == 2
        rho = spins_readdata('Density', ii, x_inds, 1, z_inds);
    end

    % set-up the figure
    all_conts = true; % do all contours exist in the domain?
    figure(18), clf
    hold on

    % loop through the contours
    for nn = 1:n_cont
        if contval(nn) < max(rho(:)) && contval(nn) > min(rho(:))
            % find contour (isopycnal)
            if isvector(gd.x)
                [cont_x, cont_y] = find_contour(gd.x(x_inds), gd.z(z_inds), rho', contval(nn));
            else
                [cont_x, cont_y] = find_contour(gd.x(x_inds,z_inds), gd.z(x_inds,z_inds), rho, contval(nn));
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
                wavelength_right(jj, nn) = 0;
                wavelength_left(jj, nn)  = 0;
                reach_end = true;
            end

            % make plot to check
            p_hand(nn) = plot(cont_x, cont_y, '.', 'Color', c_map(nn,:));
            % plot wave centre and wavelengths
            plot([1 1]*max_pos(1), [0 max_val(1)],'-', 'Color', c_map(nn,:))
            plot([1 1]*front_loc,  [0 max_val(1)],'-', 'Color', c_map(nn,:))
            plot([1 1]*back_loc,   [0 max_val(1)],'-', 'Color', c_map(nn,:))
        else
            warning(['contour ',num2str(contval(nn)),' is not in the data set, skipping.'])
            all_conts = false;
        end
    end
    % other stuff
    grid on
    title(['t=',int2str(params.plot_interval*ii),'s'])
    if all_conts
        leg = arrayfun(@(contval) ['rho = ',num2str(contval)],contval,'Uni',0);
        legend(p_hand, leg)
    end
    drawnow
    hold off

    % update reading indices
    if jj > 1
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
time = outputs*params.plot_interval;
wave_speed = 0*amplitude;
for nn = 1:n_cont
    wave_speed(:,nn) = FiniteDiff(time,1,2,false,false)*wave_center(:,nn);
end

% save data
save(filename,'time','amplitude','wave_center',...
     'wavelength_right','wavelength_left','strat_loc','wave_speed','contval');
