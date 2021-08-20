clc; clearvars; close all;
cd('../3D_101020_46')
savefnm = 'BroadSurge3D.mp4';
vid = VideoWriter(savefnm, 'MPEG-4');
vid.Quality = 100;
vid.FrameRate = 1;
open(vid);
figure(1);
clf
set(gcf, 'Units', 'centimeters'); set(gcf, 'Position', [1 1 44 20]);
set(gcf, 'PaperUnits','centimeters'); set(gcf, 'Position', [1 1 44 20]);
params = spins_params;
xlimits = [5.5 6.5];
[x, z] = spinsgrid2d;
x = squeeze(x(:, 1, :));
z = squeeze(z(:, 1, :));

xlimits_ind = find(x(:, 1)>xlimits(1) & x(:, 1)<xlimits(2));
x = x(xlimits_ind, :);
z = z(xlimits_ind, :);

for i = 85:90
    u = xz_reader('u', i, xlimits_ind, []);
    subaxis(2, 1, 1);
    pcolor(x, z, u); shading flat; cmocean('delta');
    daspect([1 1 1]);
    c1 = colorbar;
    c1.YLabel.String = 'u (m/s)';
    caxis([-.1 .1])
    rho = xz_reader('rho', i, xlimits_ind, []);
    rho = rho_converter(rho);
    subaxis(2, 1, 2)
    pcolor(x, z, rho); shading flat; cmocean('dense');
    daspect([1 1 1]);
    c = colorbar;
    c.YLabel.String = '$\rho$ (kgm$^{-3}$)';
    caxis([1026 1046]);
    
    figure_print_format(gcf);
    F = getframe(gcf);
    writeVideo(vid, F);
    
end
close(vid)
%%
clc; clearvars; close all;
x = xgrid_reader; xinds = find(x(:, 1)>5.5);

x = x(xinds,:, :);
y = ygrid_reader; y = y(xinds,:, :);
z = zgrid_reader; z = z(xinds,:, :);
zi = 59; 

disp(num2str(z(1, 1, zi)));
figure;
for i = 85:90
u = (spins_reader_new('v', i, xinds, [], []));
wh=squeeze(u(:,:,zi)); xh=squeeze(x(:,:,zi)); yh=squeeze(y(:,:,zi)); pcolor(xh,yh,wh), shading flat
cmocean('balance'); colorbar;
xlim([5.5 7])
ylim([0 .128])
pause(.2);

end
