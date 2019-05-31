% PLOT_WAVE_CHAR    Plot the data from the characterize_wave.m script
%
%  Usage: After characterize_wave has been successfully completed, run
%         this script to plot all the outputs.
%
% David Deepwell, 2017

% Set-up
load('wave_characteristics.mat')
leg = arrayfun(@(contval) ['rho = ',num2str(contval)],contval,'Uni',0);
n_cont = length(contval);
fnum = 125;

%% Create plots
% amplitude
figure(fnum)
  clf
  plot(time,amplitude*100)
  xlabel('time (s)')
  ylabel('Amplitude (cm)')
  legend(leg)
  legend('location','best')
  legend('boxoff')
  fnum = fnum+1;

% wavelengths
if n_cont == 1
figure(fnum)
  clf
  if isvector(wavelength_right)
      plot(time,[wavelength_right wavelength_left]*100)
  else
      plot(time,[wavelength_right; wavelength_left]*100)
  end
  xlabel('time (s)')
  ylabel('Wavelength (cm)')
  legend('Right','Left')
  legend('location','best')
  legend('boxoff')
  fnum = fnum+1;
else
figure(fnum)
  clf
  plot(time,wavelength_right*100)
  xlabel('time (s)')
  ylabel('Right Wavelength (cm)')
  legend(leg)
  legend('location','best')
  legend('boxoff')
  fnum = fnum+1;
figure(fnum)
  clf
  plot(time,wavelength_left*100)
  xlabel('time (s)')
  ylabel('Left Wavelength (cm)')
  legend(leg)
  legend('location','best')
  legend('boxoff')
  fnum = fnum+1;
end

% wave speed
figure(fnum)
  clf
  plot(time,wave_speed*100)
  xlabel('time (s)')
  ylabel('Wave speed (cm/s)')
  legend(leg)
  legend('location','best')
  legend('boxoff')
  fnum = fnum+1;

% wave center location
figure(fnum)
  clf
  plot(time, wave_center)
  xlabel('time (s)')
  ylabel('Wave center (m)')
  legend(leg)
  legend('location','best')
  legend('boxoff')
  fnum = fnum+1;

% background isopycnal location
figure(fnum)
  clf
  plot(time, strat_loc*100)
  xlabel('time (s)')
  ylabel('Background isopycnal location (cm)')
  legend(leg)
  legend('location','best')
  legend('boxoff')
  fnum = fnum+1;

% wavelength vs. amplitude
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

% wave speed vs. amplitude
figure(fnum)
  clf
  plot(amplitude*100, wave_speed*100,'.')
  xlabel('Amplitude (cm)')
  ylabel('Wavespeed (cm/s)')
  legend(leg)
  legend('location','best')
  legend('boxoff')
  fnum = fnum+1;

% wave area
figure(fnum)
  clf
  plot(time,amplitude.*(wavelength_left+wavelength_right)*(100^2))
  xlabel('time (s)')
  ylabel('a*(\lambda_a + \lambda_f) (cm^2)')
  legend(leg)
  legend('location','best')
  legend('boxoff')
