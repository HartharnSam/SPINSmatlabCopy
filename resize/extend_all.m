function extend_all(ii, Ly, Ny)
%EXTEND_ALL  Adjust the resolution in all dimensions of all fields.
%
%  Assumptions:
%    - spins.conf file must be present
%
%  Usage:
%    extend_all(output, spanwise_width, spanwise_gridpoints)
%
%  Inputs:
%    'output'            - output number to extend
%    'spanwise_width'    - Width of the extended dimension
%    'Ny'                - number of gridpoints in the extended dimension

%
%  Outputs:
%    n/a
%
%  Andrew Grace, 2021.

fields = find_fields(ii);

%Make new directory to store resized data
new_dir = 'extended';
mkdir(new_dir)


% read parameters
params = spins_params();
perturb = params.perturb;

Nx = params.Nx;
Nz = params.Nz;

Lx = params.Lx;
Lz = params.Lz;

for jj = 1:length(fields)
    tic
    field = fields{jj};
    fprintf('Extending field: %-8s ...',field)

    % read and expand field
    data = spins_reader(field, ii);
    
    data3D = extend(data,Ny,2);
    data3D = data3D + perturb*(2*rand(Nx,Ny,Nz) - 1);
    %Perturb data with white noise. 
    

    % write new field in new directory
    cd(new_dir)
    spins_writer([field,'.',num2str(ii)], data3D);
    cd('..')

    fprintf(' took %.4g s\n',toc)
end

%Generate a new v velocity field
v = zeros(Nx,Ny,Nz) + perturb*(2*rand(Nx,Ny,Nz) - 1);
cd(new_dir)
spins_writer(['v.',num2str(ii)], v);
cd('..')


% create grids
fprintf('Creating / writing grids ...')
tic

xg = generate_xgrid([Lx Ly Lz],[Nx Ny Nz],params.type_z);
yg = generate_ygrid([Lx Ly Lz],[Nx Ny Nz],params.type_z);
zg = generate_zgrid([Lx Ly Lz],[Nx Ny Nz],params.type_z);


% write grids in new directory
cd(new_dir)
spins_writer('xgrid', xg);
spins_writer('zgrid', zg);
spins_writer('ygrid', yg);

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

% Replace old Ly and Ny sizes and restarting information in spins.conf
cd(new_dir)
comp = computer();
restart_time = params.plot_interval*ii;
if strncmp(comp,'MAC',3)
    system(['sed -i '''' ''s/^Ly.*$/Ly = ',num2str(Ly),'/'' spins.conf']);
    system(['sed -i '''' ''s/^Ny.*$/Ny = ',num2str(Ny),'/'' spins.conf']);       
    system(['sed -i '''' ''s/^restart[[:space:]]*=.*$/restart = true/g'' spins.conf']);
    system(['sed -i '''' ''s/^restart_time.*$/restart_time = ',num2str(restart_time),'/g'' spins.conf']);
    system(['sed -i '''' ''s/^restart_sequence.*$/restart_sequence = ',num2str(ii),'/g'' spins.conf']);
    system(['sed -i '''' ''s/^restart_from_dump.*$/restart_from_dump = false/g'' spins.conf']);
else
    system(['sed -i -e ''s/^Ly.*$/Ly = ',num2str(Ly),'/g'' spins.conf']);
    system(['sed -i -e ''s/^Ny.*$/Ny = ',num2str(Ny),'/g'' spins.conf']);       
    system(['sed -i -e ''s/^restart[[:space:]]*=.*$/restart = true/g'' spins.conf']);
    system(['sed -i -e ''s/^restart_time.*$/restart_time = ',num2str(restart_time),'/g'' spins.conf']);
    system(['sed -i -e ''s/^restart_sequence.*$/restart_sequence = ',num2str(ii),'/g'' spins.conf']);
    system(['sed -i -e ''s/^restart_from_dump.*$/restart_from_dump = false/g'' spins.conf']);
end
cd('..')
fprintf('                         ... took %.4g s\n',toc)

end

