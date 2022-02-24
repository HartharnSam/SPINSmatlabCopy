clearvars;
load('all_diagnos.mat')
params = spins_params;
Max_Density_Times = interp1(all_diagnos.diagnos.Time, all_diagnos.diagnos.Max_density, [0:200]);
plot(Max_Density_Times)
hold on
yline(params.delta_rho/2)

time_inds = find(round(Max_Density_Times, 4) > params.delta_rho/2)-1;
[x z] = spinsgrid2d;
figure;
savefnm = 'V_Density_Variance.mp4';
    if ispc
        vid = VideoWriter(savefnm, 'MPEG-4');
    else
        vid = VideoWriter(savefnm);
    end
    vid.FrameRate = 1;
    open(vid);
for ii = 100
   
    rho = spins_reader_new('rho', ii, [], []);
    subplot(2, 1, 1)
    pcolor(x, z, rho); shading flat; hold on;
    caxis([-.0095 .0095])
    xlim([0 14.5])
    ylim([-.3 0])
    colorbar;
    
    subplot(2,1,  2)
    
    rho(rho<=params.delta_rho/2) = NaN;
    rho = (rho-params.delta_rho/2)/params.delta_rho;
    pcolor(x, z, rho); shading flat; hold on;
    plot(x(:, 1), z(:, 1), 'k-')
    title(num2str(ii))
    xlim([0 14.5])
    ylim([-.3 0])
    caxis([0 .2])
    colorbar;
    
    pause(.5)
    
    hold off
     figure_print_format(gcf);
        F = getframe(gcf);
        writeVideo(vid, F);
end
close(vid)