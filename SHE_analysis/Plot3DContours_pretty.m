%%
%% Countour?
clc; clearvars; close all;
savefnm = 'BroadSurge3D_Contours_pretty.mp4';
vid = VideoWriter(savefnm, 'MPEG-4');
vid.Quality = 100;
vid.FrameRate = 1;
open(vid)

%%
set(gcf, 'Units', 'centimeters'); set(gcf, 'Position', [1 1 26.5113 13.0969]);
set(gcf, 'PaperUnits','centimeters'); set(gcf, 'PaperPosition', get(gcf, 'Position'));
set(gcf, 'Color', '#262626')
%% Run the thing
for t = 85:91
    x = xgrid_reader; xinds = find(x(:, 1)>5.8 & x(:, 1)<6.5);
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
    %% Now plot it
    clf; clc;
    subaxis(1, 1, 1, 'Margin', 0);
    isos = [1030 1035 1040];
    cmp = cmocean('dense', 3);
    for j = 1:3
        
        faces = isosurface(xi, yi, zi, field_rect, isos(j));
        p(j) = patch(faces);
        zlim([-.24 0]); ylim([0 .128]); xlim([5.8 6.5])
        view(-27,  14); daspect([1 1 1]);
        axis tight
        p(j).FaceAlpha = .7;
        p(j).FaceColor = cmp(j, :);
        p(j).EdgeColor = 'none';
        hold on
    end
    hold on
    x1 = x(1, 1, 1); x2 = x(end, 1, 1);
    y1 = y(1, 1, 1); y2 = y(1, end, 1);
    z1 = z(1, 1, 1); z2 = z(end, 1, 1);
    str = '#262626';
    color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
    p1 = patch([x1 x2 x2 x1], [y1 y1 y2 y2], [z1 z2 z2 z1], color, 'EdgeColor', 'none');
    p2 = patch([x2 x1 x1 x2], [y1 y1 y1 y1], [-.3 -.3 z1 z2], color, 'EdgeColor', 'none');
    p3 = patch([x2 x2 x2 x2], [y1 y2 y2 y1], [-.3 -.3 z2 z2], color, 'EdgeColor', 'none');
    p2.FaceColor = '#262626'; p3.FaceColor = '#262626';
    p1.FaceColor = 'w'; p1.EdgeColor = '#262626'; p1.FaceAlpha = 1;
    
    % Format
    %xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
    pos = get(gca, 'Position');
    %legend(p , num2str(isos(1)), num2str(isos(2)), num2str(isos(3)), 'Location', 'west')
    %title(['t = ', num2str(t), ' s'])
    xticks([]); yticks([]); zticks([]);
    ax = gca;
    ax.XColor = '#262626';
    ax.ZColor = '#262626';
    ax.YColor = '#262626';
    figure_print_format(gcf);

    %% Write frame
    F = getframe(gcf);
    writeVideo(vid, F);
end
close(vid);
