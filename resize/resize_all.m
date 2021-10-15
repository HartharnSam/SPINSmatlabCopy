function resize_all(ii, new_grid)
%RESIZE_ALL  Adjust the resolution in all dimensions of all fields.
%
%  Assumptions:
%    - if 2D, data must be in x-z plane
%    - spins.conf file must be present
%
%  Usage:
%    resize_all(output, [new_szx new_szy new_szz])
%
%  Inputs:
%    'output'     - output number to resize
%    'new_szx'    - number of points to change to in x
%    'new_szy'    - number of points to change to in y
%    'new_szz'    - number of points to change to in z
%
%  Outputs:
%    n/a
%
%  Andrew Grace, 2021. Generalization of scripts by David Deepwell.


%Check size of new_grid to infer dimensionality of data
if length(new_grid) == 3
    %Data is 3D. Do nothing
elseif length(new_grid) == 2
    new_grid(3) = new_grid(2);
    new_grid(2) = 1;
else
    fprintf('new_grid must be of length 2 or 3');
    return;
end

% read parameters
params = spins_params();

% shorten some parameters
dt = params.plot_interval;

Nx_old = params.Nx;
Ny_old = params.Ny;
Nz_old = params.Nz;

Lx = params.Lx;
Ly = params.Ly;
Lz = params.Lz;

Nx_new = new_grid(1);
Ny_new = new_grid(2);
Nz_new = new_grid(3);

% find which fields exist at output ii
fields = find_fields(ii);

% check x sizes
mult = Nx_new/Nx_old;
if mult > 1
    ignore_x = false;
    change = 'Increasing';
    fac_x = mult;
    fprintf('%s x resolution by factor of %d\n',change,fac_x)
elseif mult < 1
    ignore_x = false;
    change = 'Decreasing';
    fac_x = 1/mult;
    fprintf('%s x resolution by factor of %d\n',change,fac_x)
else
    ignore_x = true;
    fprintf('x resolution unchanged\n')
end

% check y sizes

if Ny_old ~= 1
    mult = Ny_new/Ny_old;
    if mult > 1
        ignore_y = false;
        change = 'Increasing';
        fac_y = mult;
        fprintf('%s y resolution by factor of %d\n',change,fac_y)
    elseif mult < 1
        ignore_y = false;
        change = 'Decreasing';
        fac_y = 1/mult;
        fprintf('%s y resolution by factor of %d\n',change,fac_y)
    else
        ignore_y = true;
        fprintf('y resolution unchanged\n')
    end
else
    ignore_y = true;
end

% check z sizes
mult = Nz_new/Nz_old;
if mult > 1
    ignore_z = false;
    change = 'Increasing';
    fac_z = mult;
    fprintf('%s z resolution by factor of %d\n',change,fac_z)
elseif mult < 1
    ignore_z = false;
    change = 'Decreasing';
    fac_z = 1/mult;
    fprintf('%s z resolution by factor of %d\n',change,fac_z)
else
    ignore_z = true;
    fprintf('z resolution unchanged\n')
end

if ignore_x && ignore_y && ignore_z 
        fprintf('No resizing required.\n');
        return;
end

%Make new directory to store resized data
new_dir = 'resized';
mkdir(new_dir)

for jj = 1:length(fields)
    tic
    field = fields{jj};
    fprintf('Changing field: %-8s ...',field)

    % read and expand field
    data = spins_reader(field, ii);
    
    %Account for any unchagned sizes
    
    if ~ignore_x && ~ignore_y && ~ignore_z
        data = resize_xyz(field, data, new_grid);
    elseif ignore_x && ~ignore_y && ~ignore_z
        data = resize_y(field, data, new_grid(2));
        data = resize_z(field, data, new_grid(3));
    elseif ~ignore_x && ignore_y && ~ignore_z
        data = resize_x(field, data, new_grid(1));
        data = resize_z(field, data, new_grid(3));
    elseif ~ignore_x && ~ignore_y && ignore_z
        data = resize_x(field, data, new_grid(1));
        data = resize_y(field, data, new_grid(2));
    elseif ~ignore_x && ignore_y && ignore_z 
        data = resize_x(field, data, new_grid(1));
    elseif ignore_x && ~ignore_y && ignore_z  
        data = resize_y(field, data, new_grid(2));
    elseif ignore_x && ignore_y && ~ignore_z 
        data = resize_z(field, data, new_grid(3));
    else
        fprintf('I must have missed this case');
        return;
    end

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


if Ny_old == 1
    xg = generate_xgrid([Lx Lz],[Nx_new Nz_new],params.type_z);
    zg = generate_zgrid([Lx Lz],[Nx_new Nz_new],params.type_z);
else
    xg = generate_xgrid([Lx Ly Lz],[Nx_new Ny_new Nz_new],params.type_z);
    yg = generate_ygrid([Lx Ly Lz],[Nx_new Ny_new Nz_new],params.type_z);
    zg = generate_zgrid([Lx Ly Lz],[Nx_new Ny_new Nz_new],params.type_z);
end


% write grids in new directory
cd(new_dir)
spins_writer('xgrid', xg);
spins_writer('zgrid', zg);
if Ny_old ~=1
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

% Replace old Nx, Ny, and Nz sizes and restarting information in spins.conf
cd(new_dir)
comp = computer();
restart_time = params.plot_interval*ii;
if strncmp(comp,'MAC',3)
    system(['sed -i '''' ''s/^Nx.*$/Nx = ',num2str(Nx_new),'/'' spins.conf']);
    system(['sed -i '''' ''s/^Ny.*$/Ny = ',num2str(Ny_new),'/'' spins.conf']);   
    system(['sed -i '''' ''s/^Nz.*$/Nz = ',num2str(Nz_new),'/'' spins.conf']);    
    system(['sed -i '''' ''s/^restart[[:space:]]*=.*$/restart = true/g'' spins.conf']);
    system(['sed -i '''' ''s/^restart_time.*$/restart_time = ',num2str(restart_time),'/g'' spins.conf']);
    system(['sed -i '''' ''s/^restart_sequence.*$/restart_sequence = ',num2str(ii),'/g'' spins.conf']);
    system(['sed -i '''' ''s/^restart_from_dump.*$/restart_from_dump = false/g'' spins.conf']);
else
    system(['sed -i -e ''s/Nx\s*=\s*',num2str(Nx_old),'/Nx = ',num2str(Nx_new),'/'' spins.conf']);
    system(['sed -i -e ''s/Ny\s*=\s*',num2str(Ny_old),'/Ny = ',num2str(Ny_new),'/'' spins.conf']);  
    system(['sed -i -e ''s/Nz\s*=\s*',num2str(Nz_old),'/Nz = ',num2str(Nz_new),'/'' spins.conf']);    
    system(['sed -i -e ''s/^restart[[:space:]]*=.*$/restart = true/g'' spins.conf']);
    system(['sed -i -e ''s/^restart_time.*$/restart_time = ',num2str(restart_time),'/g'' spins.conf']);
    system(['sed -i -e ''s/^restart_sequence.*$/restart_sequence = ',num2str(ii),'/g'' spins.conf']);
    system(['sed -i -e ''s/^restart_from_dump.*$/restart_from_dump = false/g'' spins.conf']);
end
cd('..')
fprintf('                         ... took %.4g s\n',toc)

end

