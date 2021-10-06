%%
%% Countour?
clc; clearvars; close all;
savefnm = 'BroadSurge3D_Contours.mp4';
vid = VideoWriter(savefnm, 'MPEG-4');
vid.Quality = 100;
vid.FrameRate = 1;
open(vid)

%%
set(gcf, 'Units', 'centimeters'); set(gcf, 'Position', [1 1 44 20]);
set(gcf, 'PaperUnits','centimeters'); set(gcf, 'Position', [1 1 44 20]);

%x_limits = [6.2 6.8]; % For _49
x_limits = [6 6.4];
%% Run the thing
%for t = 85:91 % For _46
% for t = 95:100 % For _48
 for t = 83:90 % For _49
    x = xgrid_reader; xinds = find(x(:, 1)>x_limits(1) & x(:, 1)<x_limits(2));
    x = x(xinds,:, :);
    y = ygrid_reader; y = y(xinds,:, :);
    z = zgrid_reader; z = z(xinds,:, :);
    rho2 = rho_converter(spins_reader_new('rho', t, xinds, [], []));
    
    gd.x = squeeze(x(:, 1, :)); gd.y = y; gd.z = squeeze(z(:, 1, :));
    field_rect = nan(size(rho2));% xi = field_rect; %zi = xi;
    for i = 1:size(y, 2)
        [gd_rect, field_recttmp] = interp_onto_rect_grid(gd, squeeze(rho2(:,i, :)));
        field_recttmp(isnan(gd_rect.x)) = NaN;
        field_recttmp(isnan(gd_rect.z)) = NaN;
        field_recttmp(isnan(squeeze(y(:, 1, :)))) = NaN;
        field_rect(:, i, :) = field_recttmp;
    end
    [xi, yi, zi] = meshgrid(squeeze(gd_rect.x(:, 1)), squeeze(y(1, :, 1)), squeeze(gd_rect.z(1, :)));
    field_rect = permute(field_rect, [2, 1, 3]);
    % Now plot it
    clf; clc;
    isos = [1031 1036 1041];
    cmp = cmocean('dense', 3);
    p = gobjects(3, 1);
    for j = 1:3
        
        faces = isosurface(xi, yi, zi, field_rect, isos(j));
        p(j) = patch(faces);
        zlim([-.3 0]); ylim([0 .128]); xlim(x_limits)
        view(-28.1486,  23.4330); daspect([1 1 1]);
        axis tight
        p(j).FaceAlpha = .7;
        p(j).FaceColor = cmp(j, :);
        p(j).EdgeColor = 'none';
        hold on
    end
    hold on
    p1 = patch(x(:, :, 1), y(:, :, 1), z(:, :, 1), 'k');
    p1.FaceColor = 'k'; p1.EdgeColor = 'k'; p1.FaceAlpha = 0;
        x1 = x(1, 1, 1); x2 = x(end, 1, 1);
    y1 = y(1, 1, 1); y2 = y(1, end, 1);
    z1 = z(1, 1, 1); z2 = z(end, 1, 1);
    
    p1 = patch([x1 x2 x2 x1], [y1 y1 y2 y2], [z1 z2 z2 z1], 'k');
    %p2 = patch([x2 x1 x1 x2], [y1 y1 y1 y1], [-.3 -.3 z1 z2], 'k');
    %p3 = patch([x2 x2 x2 x2], [y1 y2 y2 y1], [-.3 -.3 z2 z2], 'k');
    % Format
    xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
    pos = get(gca, 'Position');
    legend(p , num2str(isos(1)), num2str(isos(2)), num2str(isos(3)), 'Location', 'west')
    title(['t = ', num2str(t), ' s'])
    
        figure_print_format(gcf);

    % Write frame
    F = getframe(gcf);
    writeVideo(vid, F);
end
close(vid);
