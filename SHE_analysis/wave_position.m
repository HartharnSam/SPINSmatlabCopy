%clearvars; close all; clc;

params = spins_params;
%v = VideoWriter('v3.mp4', 'MPEG-4');
%v.FrameRate = 1;
%open(v);
gd.z = zgrid_reader;
gd.x = xgrid_reader-params.L_adj;
%gdpar = spins_gridparams('Full');
%split_gdpar
% shorten parameters
%noutputs = length(dir('u.*'));
gdnames = fieldnames(gd);
params.ndims = length(gdnames);

Nx = params.Nx;
Nz = params.Nz;

startframe = 50; endframe = 80;
first_out = max(first_output(), startframe);
last_out = min(last_output(), endframe);
outputs = first_out:last_out;
noutputs = length(outputs);

strat = spins_reader_new('rho', 0, Nx/2, 1, []);
xlind = 2048;
xrind = Nx;
zbind = 1;
ztind = Nz-32;
mult = 2; % multiplier of wavelength to compute region to load
x_inds = xlind:xrind;
z_inds = zbind:ztind;

config = {'depr'};
isopyc_loc = params.pyc_loc -params.h_halfwidth + 0.01;
contval = interp1(gd.z(Nx/2,:), strat, isopyc_loc);
n_cont = length(contval);
amplitude = zeros(noutputs, n_cont);
c_map = [rgb('Black'); rgb('Red')];
time = get_output_times();
%time = (startframe:endframe)';
for jj = 1:noutputs
    ii = outputs(jj);
    
    figure(1), clf
    hold on
    % find background depth of chosen isopycnal (contval)
    
    [strat_pos, ~] = find_position(gd.z(Nx/2,:), strat, contval(1));
    
    strat_loc(jj, 1) = strat_pos;
    
    rho = spins_reader_new('rho', ii, x_inds, 1, z_inds);
    vort = spins_reader_new('vorty', ii, x_inds, 1, z_inds);
    
    %%
    if contval(1) < max(rho(:)) && contval(1) > min(rho(:))
        
        % find contour (isopycnal)
        
        [cont_x, cont_y] = find_contour(gd.x(x_inds,z_inds), gd.z(x_inds,z_inds), rho, contval(1));
        
        %cont_y = cont_y - strat_loc(jj, 1); % shift so that cont_y=0 at far field
        % if wave is a depression, then flip vertically
        
        
        % find amplitudes and locations of local maxima
        [max_val, max_pos, max_ind] = find_wave_max(cont_x, -cont_y);
        amplitude(jj, 1)   = max_val(1);
        % what to do about other components?
        
        
        % make plot to check
      
        hold on
        % plot wave centre and wavelengths
       % plot([1 1]*max_pos(1), [0 -max_val(1)]-max_val(1)/2,'-', 'Color', c_map(1,:))
        xlim([5 params.Lx])
        ylim([params.min_z 0 ])
        %plot([1 1]*front_loc,  [0 max_val(1)],'-', 'Color', c_map(1,:))
        %plot([1 1]*back_loc,   [0 max_val(1)],'-', 'Color', c_map(1,:))'
        title(outputs(jj))
        text(params.Lx-0.7, -.25, ['mid wave = ', num2str(max_pos(1)), 'm']);
        plot([params.Lx, params.Lx-params.hill_height/params.hill_slope]-params.L_adj, params.min_z+[params.hill_height 0], 'k-', 'linewidth', .2);
        c = contour(gd.x(x_inds,z_inds), gd.z(x_inds,z_inds), vort, [-5 -4 4 5]);
        caxis([-6 6]);
        colormap(brighten(darkjet, .8));
        %alpha(c, .5)
        colorbar
        
        p_hand(1) = plot(cont_x, cont_y, '.', 'Color', c_map(1,:));
        %F = getframe(gcf);
        %writeVideo(v, F);
    else
        warning(['contour ',num2str(contval(1)),' is not in the data set, skipping.'])
        all_conts = false;
    end
%     inds = max_ind(1):length(cont_x);
%     if length(inds) <= 1
%         disp(['Wave has reached the tank end, skipping output '...
%             ,num2str(ii+1),' and all after'])
%         completion(jj, noutputs)
%         break
%     end
    
    %(2); clf;
    %plot([params.Lx, params.Lx-params.hill_height/params.hill_slope]-params.L_adj, params.min_z+[params.hill_height 0], 'k-', 'linewidth', .2);
    %hold on
    
end
%close(v)