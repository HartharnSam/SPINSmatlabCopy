function [] = make_movie_bare()
% MAKE_MOVIE_BARE     Basic movie making script.
%   Uses ffmpeg to stitch frames together, so a working copy 
%   is necessary.
%
%  Usage:
%    make_movie_bare()
%
%  Inputs:
%   none
%
%  Outputs:
%    none
%
%  David Deepwell, 2017

% video frame rate
framerate = 3;  % frame rate lower than 3 does not work with quicktime

% set-up directory
if exist('tmp_figs','dir')
    rmdir('tmp_figs','s') % wipe directory clear before populating
end
mkdir('tmp_figs')

% outputs to use in movie (default here is all of them)
first_out = first_output();
last_out  =  last_output();
outputs   = first_out:last_out;

% do loop
for ii = outputs
    % make plot (put function for making a frame here)

    % save output figure
    cd('tmp_figs')
    filename = ['tmp_',num2str(ii,'%03d')];
    print_figure(filename, 'format', 'png')
    cd('..')

    completion(ii-first_out+1, length(outputs))
end

% create directory to store movies
if ~(exist('movies','dir')==7)
    mkdir movies
end
% make movie and remove temporary files and folder
cd('tmp_figs')
filename = 'movie';
status = system(['ffmpeg -r ',num2str(framerate),...
         ' -start_number ',num2str(first_out),' -i tmp_%03d.png',...
         ' -r ',num2str(framerate),' -y -pix_fmt yuv420p -q 1',...
         ' -vf scale=-1:600 ../movies/',filename,'.mp4']);
if status ~= 0
    disp([filename,'.mp4 was possibly not redered correctly.'])
end
cd('..')
rmdir('tmp_figs','s')
