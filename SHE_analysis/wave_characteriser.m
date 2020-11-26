clc ; clearvars; close all;
params = spins_params;
%end_point = params.Lx * 2/3;
%start_point = params.Lx * 1/4;
threshold = 1/.007;
%x_points = linspace(start_point, end_point, nt);
times = [25:50];

for i = 1:length(times)
    %[time, z, data] = spins_hovmoller('rho',  x_points(i));
    ii = times(i);
    xi = xgrid_reader; 
    z = zgrid_reader;
    data = spins_reader_new('rho', ii);
    M = contour(xi, z, data, [0 0], 'k-');
    M = M(:, M(1, :)<params.Lx - params.hill_height/params.hill_slope);
    
    [t, ind] = sort(M(1, :));
    t = flip(t);
    x = M(2, ind);
    x = flip(x);
    x = smoothdata(x, 'movmedian', 5);
    % %find crossing points
    Mn = round((x +.15) * threshold);
    %plot(M(1, 2:end), Mn)
    cross_points = diff(sign(Mn)); %A matrix where a value of 2 represents a -ve to +ve change between two points
    point_1 = find(Mn~=0, 1, 'first');
    hold on
    plot(t(point_1), x(point_1), 'kx'); % Plot where first wave starts
    
    %%
    indx_up = find(cross_points >=1, 1, 'first'); %Indexes points where a negative->positive change occurs
    indx_down = find(cross_points <=-1, 1, 'first');
    plot(t(indx_up), x(indx_up), 'kx'); % Plot where first wave ends
    plot(t(indx_down), x(indx_down), 'kx'); % Plot where first wave ends
    
    t_up = t(indx_up);
    t_dwn = t(indx_down);
    
    %%
    [amp(i), amp_ind(i)] = max(abs(x(indx_down:indx_up)));
    amp(i) = amp(i) - .15;
    amp_ind(i) = amp_ind(i) + indx_down;
    amp_time(i) = t(amp_ind(i));
    plot(t(amp_ind(i)), -amp(i)-.15, 'rx')
    
    %% Second find 0 crossing
    cross_points = diff(sign((x(indx_down:indx_up) +.15 +amp(i)/2)));
    %A matrix where a value of 2 represents a -ve to +ve change between two points
    front_indx_point = indx_down+find(cross_points <=-1, 1, 'first');
    front_wl_point = t(front_indx_point);
    rear_indx_point = indx_down+find(cross_points >= 1, 1, 'first');
    rear_wl_point(i) = t(rear_indx_point);
    plot(rear_wl_point(i), x(rear_indx_point), 'bx');
    plot(front_wl_point, x(front_indx_point), 'bx');
    
    right_wl_dist(i) = abs(rear_wl_point(i) - amp_time(i));
    wl_time_dist(i) = abs(rear_wl_point(i) - front_wl_point);
    disp([num2str(i), ' of ', num2str(length(times))])
    %set(gca, 'XDir', 'normal');
    xlabel('x')
    ylabel('z')
    if ii == times(1) || ii == times(end)
    print(['C:\Users\samha\OneDrive - Newcastle University\Project\Shoal_Core\figures\Continuous_Contour_', num2str(i), '.png'], '-dpng');
    end
    hold off
end

c = diff(amp_time)./diff(times);
wavelength = wl_time_dist(2:end);
right_full_wavelength = right_wl_dist*2;
mean_wavelength = mean(wavelength);
mean_right_full_wavelength = mean(right_full_wavelength);
WaveStats.meanAmp = mean(amp);
WaveStats.meanWaveSpeed = mean(c);
WaveStats.meanWavelength = mean_wavelength;
save('wavestats.mat', 'WaveStats');