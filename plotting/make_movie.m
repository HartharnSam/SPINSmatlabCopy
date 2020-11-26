function [] = make_movie(varname)
% MAKE_MOVIE     Make a movie tracking the wave in the current directory.
%   make_movie uses ffmpeg to stitch frames together, so a working copy 
%   of that is necessary.
%
%  Usage:
%    make_movie('rho')
%
%  Inputs:
%    'varname'   - the variable to make a movie of.
%                  Allowed names - Anything readable by spins_plot2d:
%                      rho, u, KE, Mean KE, SD KE, ...
%
%  Outputs:
%    none
%
%  David Deepwell, 2018

% get grid and parameters
params = spins_params;

% colorbar option
if strcmp(varname,'rho')
    clim = [-1 1]*params.delta_rho/2;
elseif strncmp(varname,'dye',3) || strcmp(varname,'tracer')
    clim = [0 1];
    if strcmp(varname,'dye1b')
        clim = [-1 0];
    end
else
    clim = 'auto';
end
% axis options
zm = params.Lz/2 + params.min_z;
zB = zm - 0.07;
zT = zm + 0.07;
xL = -0.4;
xR = 0.1;
width = xR - xL;

% video frame rate
framerate = 3;  % frame rate lower than 3 does not work with quicktime

% set-up directory and figure
if exist('tmp_figs','dir')
    rmdir('tmp_figs','s') % wipe directory clear before populating
end
mkdir('tmp_figs')

% outputs to use in movie (default here is all of them)
first_out = first_output(varname);
last_out = last_output(varname);
outputs = first_out:last_out;

% do loop
for ii = outputs
    % find region
    ax = wave_region(ii,'type','fixed','x',[xL xR],'z',[zB zT]);
    if ax(2)-ax(1) < width
        if ax(2) > params.Lx
            ax(1) = ax(2) - width;
        elseif ax(1) < 0
            ax(2) = ax(1) + width;
        end
    end

    % make plot
    if strcmp(clim, 'auto')
        spins_plot2d(varname,ii,'axis',ax);
    else
        spins_plot2d(varname,ii,'axis',ax,'clim',clim,'trim',true);
    end
    if strncmp(varname,'dye',3)
        colormap(cmocean('tempo'))
        if strcmp(varname,'dye1b')
            colormap(cmocean('-tempo'))
        end
    end
    % adjust to defaults
    figure_defaults()

    % save output figure
    cd('tmp_figs')
    filename = ['tmp_',num2str(ii,'%03d')];
    print_figure(filename, 'format', 'png','size',[8 3])
    cd('..')

    completion(ii-first_out+1, length(outputs))
end

% create directory to store movies
if ~(exist('movies','dir')==7)
    mkdir movies
end
% make movie and remove temporary files and folder
cd('tmp_figs')
filename = strrep(varname, ' ','_');
status = system(['ffmpeg -r ',num2str(framerate),...
         ' -start_number ',num2str(first_out),' -i tmp_%03d.png',...
         ' -r ',num2str(framerate),' -y -pix_fmt yuv420p -q 1',...
         ' -vf scale=-1:600 ../movies/',filename,'.mp4']);
if status ~= 0
    disp([varname,'.mp4 was possibly not rendered correctly.'])
end
cd('..')
rmdir('tmp_figs','s')
