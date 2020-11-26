% PLOT_WAVE_CHAR    Plot the data from the characterize_wave.m script
%
%  Usage: After characterize_wave has been successfully completed, run
%         this script to plot all the outputs.
%
% David Deepwell, 2017

% Set-up
if ~isfile('wave_characteristics.mat')
    characterize_wave;
end
load('wave_characteristics.mat')
load('wavestats.mat', 'WaveStats'); endslope = WaveStats.endSlope;
if isempty(endslope)
    endslope = length(time)
end

leg = arrayfun(@(contval) ['rho = ',num2str(contval)],contval,'Uni',0);
n_cont = length(contval);
fnum = 125;
n_plots = 7;
%% Create plots
%% amplitude
figure(fnum)
set(gcf, 'DefaultAxesXDir', 'normal');
set(gcf, 'Units', 'centimeters'); set(gcf, 'Position', [1 1 20 20]);
set(gcf, 'PaperUnits','centimeters'); set(gcf, 'Position', [1 1 20 20]);

clf
subplot(n_plots, 1, 1)
plot(time,amplitude*100)
hold on
plot([time(endslope) time(endslope)], [0 max(amplitude)*100], 'k-')
%xlabel('time (s)')
xticks([]);
ylabel({'Amplitude', '(cm)'})
legend(leg)
legend('location','best')
legend('boxoff')
fnum = fnum+1;

%% wavelengths
if n_cont == 1
    %figure(fnum)
    %  clf
    subplot(n_plots, 1, 2);
    if isvector(wavelength_right)
        plot(time,[wavelength_right wavelength_left]*100)
    else
        plot(time,[wavelength_right; wavelength_left]*100)
    end
    hold on
    plot([time(endslope) time(endslope)], [0 max(wavelength_right)*100], 'k-')
    %xlabel('time (s)')
    xticks([]);
    ylabel({'Wavelength', '(cm)'})
    legend('Right','Left')
    legend('location','best')
    legend('boxoff')
    fnum = fnum+1;
    
else
    %figure(fnum)
    %  clf
    subplot(n_plots, 1, 2)
    plot(time,wavelength_right*100)
       hold on
    plot([time(endslope) time(endslope)], [0 max(wavelength_right)*100], 'k-')

    xlabel('time (s)')
    %ylabel('Right Wavelength (cm)')
    legend(leg)
    legend('location','best')
    legend('boxoff')
    fnum = fnum+1;
    %figure(fnum)
    %clf
    
    plot(time,wavelength_left*100)
    xlabel('time (s)')
    ylabel('Wavelength (cm)')
    legend(leg)
    legend('location','best')
    legend('boxoff')
    fnum = fnum+1;
end

%% Total wavelengths
if n_cont == 1
    %figure(fnum)
    %  clf
    subplot(n_plots, 1, 3);
    if isvector(wavelength_right)
        plot(time,[wavelength_right+wavelength_left]*100)
    else
        plot(time,[wavelength_right+wavelength_left]*100)
    end
    hold on
    plot([time(endslope) time(endslope)], [0 max(wavelength_right+wavelength_left)*100], 'k-')
    %xlabel('time (s)')
    xticks([]);
    ylabel({'Wavelength', '(cm)'})
    legend('Right','Left')
    legend('location','best')
    legend('boxoff')
    fnum = fnum+1;
    
else
    %figure(fnum)
    %  clf
    subplot(n_plots, 1, 3)
    plot(time,wavelength_right*100)
       hold on
    plot([time(endslope) time(endslope)], [0 max(wavelength_right)*100], 'k-')

    xlabel('time (s)')
    %ylabel('Right Wavelength (cm)')
    legend(leg)
    legend('location','best')
    legend('boxoff')
    fnum = fnum+1;
    %figure(fnum)
    %clf
    
    plot(time,wavelength_left*100)
    xlabel('time (s)')
    ylabel('Wavelength (cm)')
    legend(leg)
    legend('location','best')
    legend('boxoff')
    fnum = fnum+1;
end

%% wave speed
%figure(fnum)
%clf
subplot(n_plots, 1, 4)
plot(time,wave_speed*100)
hold on
plot([time(endslope) time(endslope)], [-1 1]*max(wave_speed)*100, 'k-')

xlabel('time (s)')
ylim([-1 1]*max(wave_speed)*100);
ylabel({'Wave speed', '(cm/s)'})
legend(leg)
legend('location','best')
legend('boxoff')
fnum = fnum+1;

%% wave center location
% figure(fnum)
% clf
subplot(n_plots, 1, 5)
plot(time, wave_center)
hold on
plot([time(endslope) time(endslope)], [0 max(wave_center)], 'k-')
%xlabel('time (s)')
xticks([]);
ylabel({'Wave center', '(m)'})
legend(leg)
legend('location','best')
legend('boxoff')
fnum = fnum+1;

%% background isopycnal location
subplot(n_plots, 1, 6)
plot(time, strat_loc*100)
xlabel('time (s)')
ylabel({'Background', 'isopycnal', 'location', '(cm)'})
legend(leg)
legend('location','best')
legend('boxoff')
fnum = fnum+1;

%% wave area
subplot(n_plots, 1, 7)
plot(time,amplitude.*(wavelength_left+wavelength_right)*(100^2))
xlabel('time (s)')
ylabel({'a*(\lambda_a + \lambda_f)', '(cm^2)'})
legend(leg)
legend('location','best')
legend('boxoff')

print('timeWaveStats.png', '-dpng');
%% wavelength vs. amplitude
figure(fnum)
clf
plot(amplitude*100, wavelength_right*100,'.',...
    amplitude*100, wavelength_left*100, '.')
xlabel('Amplitude (cm)')
ylabel('Wavelength (cm)')
legend('Right','Left') % if more than 1 contour used this needs to be adjusted
legend('location','best')
legend('boxoff')
fnum = fnum+1;
print('WavelengthVsAmp.png', '-dpng');

%% wave speed vs. amplitude
figure(fnum)
clf
plot(amplitude*100, wave_speed*100,'.')
xlabel('Amplitude (cm)')
ylabel('Wavespeed (cm/s)')
legend(leg)
legend('location','best')
legend('boxoff')
fnum = fnum+1;
print('SpeedVsAmp.png', '-dpng');


