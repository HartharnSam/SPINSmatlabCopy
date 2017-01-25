function [] = make_movie_bare()
% MAKE_MOVIE_BARE     Make a movie tracking the wave in the current directory.
%   make_movie uses ffmpeg to stitch frames together, so a working copy 
%   of that is necessary.
%
%  Usage:
%    make_movie()
%
%  Inputs:
%   none
%
%  Outputs:
%    none

% video frame rate
framerate = 3;  % frame rate lower than 3 does not work with quicktime

% set-up directory and figure
if exist('tmp_figs','dir')
    rmdir('tmp_figs','s') % wipe directory clear before populating
end
fig_hand = figure(1);
set(fig_hand, 'Units', 'Inches')
pos = get(fig_hand, 'Position');
set(fig_hand, 'PaperPositionMode','Auto',...
              'PaperUnits', 'Inches',...
              'PaperSize', [pos(3) pos(4)]);

% outputs to use in movie (default here is all of them)
gdpar = spins_gridparams('Vector',false); % get parameters
noutputs = gdpar.params.noutputs;
first_out = first_output();
last_out = first_out + noutputs - 1;
outputs = first_out:last_out;

% do loop
for ii = outputs
    % make plot (put function for making a frame here)

    % save output figure
    cd('tmp_figs')
    filename = ['tmp_',num2str(ii,'%03d')];
    print(fig_hand, filename, '-dpng', '-r500')
    cd('..')

    completion(ii-first_out+1, length(outputs))
end

% create directory to store movies
if ~(exist('movies','dir')==7)
    mkdir movies
end
% make movie and remove temporary files and folder
cd('tmp_figs')
filename = strrep(var, ' ','_');
status = system(['ffmpeg -r ',num2str(framerate),...
         ' -start_number ',num2str(first_out),' -i tmp_%03d.png',...
         ' -r ',num2str(framerate),' -y -pix_fmt yuv420p -q 1',...
         ' -vf scale=-1:600 ../movies/',filename,'.mp4']);
if status ~= 0
    disp([var,'.mp4 was possibly not redered correctly.'])
end
cd('..')
rmdir('tmp_figs','s')
