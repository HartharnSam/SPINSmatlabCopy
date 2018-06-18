function resize_all_y(ii, Ny_new)
% RESIZE_ALL_Y  Adjust the resolution in y of all fields.
%
%  Assumptions:
%    - spins.conf file must be present
%
%  Usage:
%    resize_all_y(5, 512)
%
%  Inputs:
%    'ii'        - output number to resize
%    'Ny_new'    - number of points to change to
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
Ny_old = params.Ny;
Nz = params.Nz;
Lx = params.Lx;
Ly = params.Ly;
Lz = params.Lz;
dx = Lx/Nx;
dz = Lz/Nz;

% check sizes
mult = Ny_new/Ny_old;
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

fprintf('%s y resolution by factor of %d\n',change,fac)
% loop over fields
for jj = 1:length(fields)
    tic
    field = fields{jj};
    fprintf('Changing field: %-8s ...',field)

    % read and expand field
    data = spins_reader(field, ii);
    data = resize_y(field, data, Ny_new);

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
dy_new = Ly/Ny_new;
x1d = dx/2:dx:Lx;
z1d = dz/2:dz:Lz;
y1d = dy_new/2:dy_new:Ly;
xg = bsxfun(@times, ones(Nx, Ny_new, Nz), x1d');
yg = bsxfun(@times, ones(Nx, Ny_new, Nz), y1d);
zg = bsxfun(@times, ones(Nx, Ny_new, Nz), reshape(z1d,1,1,Nz));

% write grids in new directory
cd(new_dir)
spins_writer('xgrid', xg);
spins_writer('ygrid', yg);
spins_writer('zgrid', zg);
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

% Replace old Ny size and restarting information in spins.conf
cd(new_dir)
comp = computer();
restart_time = params.plot_interval*ii;
if strncmp(comp,'MAC',3)
    system(['sed -i '''' ''s/^Ny.*$/Ny = ',num2str(Ny_new),'/'' spins.conf']);
    system(['sed -i '''' ''s/^restart[[:space:]]*=.*$/restart = true/g'' spins.conf']);
    system(['sed -i '''' ''s/^restart_time.*$/restart_time = ',num2str(restart_time),'/g'' spins.conf']);
    system(['sed -i '''' ''s/^restart_sequence.*$/restart_sequence = ',num2str(ii),'/g'' spins.conf']);
    system(['sed -i '''' ''s/^restart_from_dump.*$/restart_from_dump = false/g'' spins.conf']);
else
    system(['sed -i -e ''s/Ny\s*=\s*',num2str(Ny_old),'/Ny = ',num2str(Ny_new),'/'' spins.conf']);
    system(['sed -i -e ''s/^restart[[:space:]]*=.*$/restart = true/g'' spins.conf']);
    system(['sed -i -e ''s/^restart_time.*$/restart_time = ',num2str(restart_time),'/g'' spins.conf']);
    system(['sed -i -e ''s/^restart_sequence.*$/restart_sequence = ',num2str(ii),'/g'' spins.conf']);
    system(['sed -i -e ''s/^restart_from_dump.*$/restart_from_dump = false/g'' spins.conf']);
end
cd('..')
fprintf('                         ... took %.4g s\n',toc)
