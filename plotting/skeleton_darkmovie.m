clc; clearvars; close all

[x, z] = spinsgrid2d;

vid = VideoWriter('../test_titlevideo.mp4', 'MPEG-4');
vid.Quality = 100; % avoids aggressive compression
vid.FrameRate = 4;
open(vid);


clrmap = linspace(255, 23, 255);
clrmap = repmat(clrmap, 3, 1)'/255;
axes('Position', [0 0 1 1]);

set(gcf, 'Units', 'centimeters', 'Position', [0 0 33 13]);
for ii = 580:680
    rho = spins_reader_new('rho', ii);

    pcolor(x, z, rho);
    ylim([.23 .3])
    colormap(clrmap);
    drawnow; pause(0.1); % This seems to be needed to make sure it's done it's thing before "getting frame"
    F = getframe(gcf);
    writeVideo(vid, F);

end
close(vid);
