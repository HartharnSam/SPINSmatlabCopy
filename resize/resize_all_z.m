function resize_all_z(ii, Nz_new)
% RESIZE_ALL_Z  Adjust the resolution in z of all fields.
%
%  Assumptions:
%    - if 2D, data must be in x-z plane
%    - spins.conf file must be present
%
%  Usage:
%    resize_all_z(5, 512)
%
%  Inputs:
%    'ii'        - output number to resize
%    'Nz_new'    - number of points to change to
%
%  Outputs:
%    n/a
%
%  David Deepwell, 2018

% read parameters
params = spins_params();

% shorten some parameters
dt = params.plot_interval;
Nx = params.Nx;
Nz_old = params.Nz;
Lx = params.Lx;
Lz = params.Lz;
dx = Lx/Nx;
try
    Ny = params.Ny;
    Ly = params.Ly;
    dy = Ly/Ny;
catch
    Ny = 1;
end

% check sizes
mult = Nz_new/Nz_old;
if mult > 1
    change = 'Increasing';
    fac = mult;
else
    change = 'Decreasing';
    fac = 1/mult;
end

% find which fields exist at output ii
fields = find_fields(ii);

% make directory to store resized data
new_dir = 'resized';
mkdir(new_dir)

fprintf('%s z resolution by factor of %d\n',change,fac)
% loop over fields
for jj = 1:length(fields)
    tic
    field = fields{jj};
    fprintf('Changing field: %-8s ...',field)

    % read and expand field
    data = spins_reader(field, ii);
    data = resize_z(field, data, Nz_new);

    % write new field in new directory
    cd(new_dir)
    spins_writer([field,'.',num2str(ii)], data);
    cd('..')

    fprintf(' took %.4g s\n',toc)
end
clear data

% create grids
fprintf('Creating / writing grids ...')
tic
dz_new = Lz/Nz_new;
x1d = dx/2:dx:Lx;
z1d = dz_new/2:dz_new:Lz;
if Ny ~= 1
    y1d = dy/2:dy:Ly;
    xg = bsxfun(@times, ones(Nx, Ny, Nz_new), x1d');
    yg = bsxfun(@times, ones(Nx, Ny, Nz_new), y1d);
    zg = bsxfun(@times, ones(Nx, Ny, Nz_new), reshape(z1d,1,1,Nz_new));
else
    xg = bsxfun(@times, ones(Nx, Nz_new), x1d');
    zg = bsxfun(@times, ones(Nx, Nz_new), z1d);
end

% write grids in new directory
cd(new_dir)
spins_writer('xgrid', xg);
spins_writer('zgrid', zg);
if Ny ~= 1
    spins_writer('ygrid', yg);
end
cd('..')
fprintf(' took %.4g s\n',toc)

% move spins.conf, case file, executable, and submit script
fprintf('Copying case files and updating spins.conf\n')
tic
system(['cp spins.conf ',new_dir]);
system(['cp *.cpp ',new_dir]);
system(['cp *.x ',new_dir]);
try
    system(['cp submit.sh ',new_dir]);
end

% Replace old Nz size and restarting information in spins.conf
cd(new_dir)
comp = computer();
restart_time = params.plot_interval*ii;
if strncmp(comp,'MAC',3)
    system(['sed -i '''' ''s/^Nz.*$/Nz = ',num2str(Nz_new),'/'' spins.conf']);
    system(['sed -i '''' ''s/^restart[[:space:]]*=.*$/restart = true/g'' spins.conf']);
    system(['sed -i '''' ''s/^restart_time.*$/restart_time = ',num2str(restart_time),'/g'' spins.conf']);
    system(['sed -i '''' ''s/^restart_sequence.*$/restart_sequence = ',num2str(ii),'/g'' spins.conf']);
    system(['sed -i '''' ''s/^restart_from_dump.*$/restart_from_dump = false/g'' spins.conf']);
else
    system(['sed -i -e ''s/Nz\s*=\s*',num2str(Nz_old),'/Nz = ',num2str(Nz_new),'/'' spins.conf']);
    system(['sed -i -e ''s/^restart[[:space:]]*=.*$/restart = true/g'' spins.conf']);
    system(['sed -i -e ''s/^restart_time.*$/restart_time = ',num2str(restart_time),'/g'' spins.conf']);
    system(['sed -i -e ''s/^restart_sequence.*$/restart_sequence = ',num2str(ii),'/g'' spins.conf']);
    system(['sed -i -e ''s/^restart_from_dump.*$/restart_from_dump = false/g'' spins.conf']);
end
cd('..')
fprintf('                         ... took %.4g s\n',toc)
