% Skeleton movie plotting script
% Basic script to efficiently produce SPINS movies in Matlab
clc; clearvars; close all; 

fig = gcf; 
fig.Position = groot().MonitorPositions(end, :); % Set figure to full screen

% Initialise video writer
vid = VideoWriter('FILENAME.mp4', 'MPEG-4');
vid.Quality = 100; % avoids aggressive compression 
vid.FrameRate = 4;
open(vid);

t1 = 0; % start time
t2 = 50; %end time

% Single Panel
ax = axes;

% Read in SPINS griddata
[x, z] = spinsgrid2d;
params = spins_params;

for ii = t1:t2
    % Run some scripts that make a plot
    rho = spins_reader_new('rho', ii);
    pcolor(x, z, rho);
    
    if ii == t1 % Set all formatting on the first iteration ONLY 
        %- this massively increases efficiency
        
        %caxis(ax(1), [-1 1].*params.delta_rho);
        colormap(ax(1), cmocean('dense'));
        axis(ax(1), 'tight');
        xlabel(ax(1), "$x (m)$")
        ylabel(ax(1), "$z (m)$");
        xlim(ax(1),[0 params.Lx]+params.min_x);
    end
    set(ax(1), 'NextPlot', 'replacechildren')
    
    drawnow; pause(0.1); % This seems to be needed to make sure it's done it's thing before "getting frame"
    figure_print_format(fig);
    F = getframe(fig);
    writeVideo(vid, F);
    completion(ii-t1, t2-t1);
end
close(vid);
