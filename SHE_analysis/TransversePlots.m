%% Pcolor of plot at z coordinates
clc; clearvars; close all;
xlimits = [6.0 6.5];
x = xgrid_reader; xinds = find(x(:, 1)> xlimits(1) & x(:, 1)<xlimits(2));
t = 91; % _46
t = 98; % _48
t = 89; % _49

x = x(xinds,:, :);
y = ygrid_reader; y = y(xinds,:, :);
z = zgrid_reader; z = z(xinds,:, :);
u = (spins_reader_new('u', t, xinds, [], []));
zi = 30;
subaxis(3, 1, 1, 'Holdaxis')
title(['t = ', num2str(t), 's'])

pcolor(squeeze(x(:, 1, :)), squeeze(z(:, 1, :)), squeeze(u(:, 1, :))); shading flat; cmocean('delta');
hold on; plot(x(:, 1, zi), z(:, 1, zi),'k--', 'LineWidth', 0.5)
plot(x(:, 1, 1), z(:, 1, 1), 'k-', 'LineWidth', .5)
xlim([xlimits])
ylim([-.2 0]); 
ylabel('z [m]');
%caxis([1026 1046]);
caxis([-.1 .1])
c = colorbar; c.YLabel.String = '$u$';
daspect([1 1 1]);
set(gca, 'Position', [.1 .4597 .7504 .5449])

subaxis(3, 1, 2, 'Holdaxis')
xh=squeeze(x(:,:,1)); yh=squeeze(y(:,:,1)); 
wh=squeeze(u(:,:,zi)); 
pcolor(xh, yh, wh); shading flat
%surf(xh,yh,zh, wh, 'FaceAlpha', .5), shading flat
%daspect([1 .08 1])
ylabel('y [m]'); zlabel('z');
xlim([xlimits])
ylim([0 .128]); zlim([-.3 0])
%caxis([1026 1046]);
caxis([-.1 .1])
hold on
zh = squeeze(z(:, :, 1));
surf(xh, yh, zh, 'FaceColor', 'k');
cmocean('delta')
c = colorbar; c.YLabel.String = '$u$';
daspect([1 1 1]);
set(gca, 'Position', [.1 .2785 .755 .2333]);

subaxis(3, 1, 3, 'Holdaxis')
rho = rho_converter(spins_reader_new('rho', t, xinds, [], zi));
pcolor(xh, yh, rho); shading flat
c = colorbar; c.YLabel.String = '$\rho$';
cmocean('dense')
xlim([xlimits])
ylim([0 .128]); zlim([-.3 0])
daspect([1 1 1]);
set(gca, 'Position', [.1 .0346 .7552 .2333]);

figure_print_format(gcf);
print(['3Panel3DView_', num2str(t), 's.png'], '-dpng');
