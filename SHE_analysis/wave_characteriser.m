function WaveStats = wave_characteriser(time_window)
%FUNCTION_NAME - track the wave core assuming it moves from left to right
%For use with  continuous stratification only
%
% Inputs:
%    time_window - [OPTIONAL 1x2 vector] - Frame numbers to characterize wave based on (e.g. [5
%    50]). Default [0 last_output]
%
% Outputs:
%    WaveStats - Structure containing key wave measurements (amplitude,
%    wavelength, speed, time it hits the slope)
%
% Other m-files required: spins_params, spins_reader_new, FiniteDiff
% Subfunctions: none
% MAT-files required: none
%
% See also: CHARACTERIZE_WAVE
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 02-Dec-2020; Last revision: 02-Dec-2020
% MATLAB Version: 9.9.0.1467703 (R2020b)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------clc ; clearvars; close all;
% read grid and parameters
params = spins_params;
Nx = params.Nx;
xi = xgrid_reader;
zi = zgrid_reader;

threshold = 1/.007;
if nargin<1
    time = 25:50;
else
    time = time_window(1):time_window(2);
end

%% Set up empty vectors
noutputs = length(time);
amplitude = zeros(noutputs, 1);
wave_center = amplitude;
wavelength_right = amplitude;
wavelength_left = amplitude;
strat_loc = amplitude;
reach_end = false;

%% Find background stratification
strat = spins_reader_new('rho', 0, Nx/2, 1, []);
isopyc_loc = params.pyc_loc;
contval = interp1(zi(Nx/2, :), strat, isopyc_loc);

for i = 1:length(time)
    if reach_end
        disp(['Wave has reached the tank end, skipping output '...
            ,num2str(i+1),' and all after'])
        continue
    end
    ii = time(i);
    % Load in grids and data
    data = spins_reader_new('rho', ii);
    
    % Find contour, and crop data to just the x region of interest
    M = contour(xi, zi, data, [contval contval], 'k-');
    M = M(:, M(1, :)<params.Lx - params.hill_height/params.hill_slope);
    
    % Put the contour data in order
    [x_cont, ind] = sort(M(1, :));
    x_cont = flip(x_cont);
    z_cont = M(2, ind);
    z_cont = flip(z_cont);
    
    % Smooth the contour to remove small numerical-linked wobbles on the
    % line
    z_cont = smoothdata(z_cont, 'movmean', 50);
    
    % %find crossing points
    Mn = round((z_cont - params.pyc_loc) * threshold); % displaces the contour to be 0 at the pyc_loc, and then puts a threshold so that small variation ignored
    
    cross_points = diff(sign(Mn)); %A matrix where a value of 2 represents a -ve to +ve change between two points
    
    indx_up = find(cross_points >=1); %Indexes points where a negative->positive change occurs
    indx_down = find(cross_points <=-1);
    x_up = x_cont(indx_up);
    x_dwn = x_cont(indx_down);
    
    hold on
    plot(x_up, z_cont(indx_up), 'bx'); % Plot where first wave ends
    plot(x_dwn, z_cont(indx_down), 'rx'); % Plot where first wave ends
    
    %%
    [amplitude(i), max_ind] = max(abs(z_cont(indx_down(1):indx_up(1))));
    amplitude(i) = amplitude(i) + params.pyc_loc; % Gives the amplitude of the wave
    max_ind = max_ind + indx_down(1);
    wave_center(i) = x_cont(max_ind);
    plot(wave_center(i), -amplitude(i)-.15, 'yx')
    
    %% Second find half_amp crossing
    zero_cross_line = z_cont(indx_down(1):indx_up(1)) -params.pyc_loc +amplitude(i)/2;
    cross_points = diff(sign(zero_cross_line));
    
    %A matrix where a value of 2 represents a -ve to +ve change between two points
    
    % Identify front and back halfwavelength points for first wave
    front_halfamp_indx = indx_down(1) + find(cross_points <= -1, 1, 'first');
    front_halfamp_pos = x_cont(front_halfamp_indx);
    wavelength_right(i) = - diff(x_cont([front_halfamp_indx max_ind]));
    
    back_halfamp_indx = indx_down(1) + find(cross_points >= 1, 1, 'first');
    back_halfamp_pos = x_cont(back_halfamp_indx);
    wavelength_left(i) = - diff(x_cont([max_ind back_halfamp_indx]));
    
    plot(front_halfamp_pos, z_cont(front_halfamp_indx), 'r*');
    plot(back_halfamp_pos, z_cont(back_halfamp_indx), 'b*');
    
    disp([num2str(i), ' of ', num2str(length(time))])
    xlabel('x')
    ylabel('z')
    
    hold off
    [strat_pos, ~] = find_position(zi(Nx/2,:), strat, contval);
    strat_loc(i) = strat_pos;
    
    if params.Lx - wave_center(i) < (wavelength_left(i)+wavelength_right(i))*1.15
        reach_end = true;
    end
    drawnow
end

wave_speed = FiniteDiff(time,1,2,false,false)*wave_center;

end_of_slope = params.Lx-((params.hill_height/params.hill_slope)+params.hill_end_dist);
read_inds = find(wave_center>(params.L_adj*1.15) & wave_center<end_of_slope);

WaveStats.endSlope = time(max(read_inds));

WaveStats.meanAmp = mean(amplitude);
WaveStats.meanWavelength = mean(wavelength_right + wavelength_left)/2;
WaveStats.meanWaveSpeed = mean(wave_speed);

save('wave_characteristics.mat', 'contval', 'isopyc_loc', 'amplitude',...
    'wave_center', 'wavelength_right', 'wavelength_left', 'time', 'strat_loc',...
    'wave_speed', 'WaveStats');