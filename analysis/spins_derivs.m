function data = spins_derivs(derivative, ii, save)
% Calculates derivatives offline. Currently only set up for 2D
% calculate dx dy dz as appropriate
% switch for derivatives (dissipation, vorticity, 
%
% Other m-files required: spins_params; spins_reader_new; get_grad2,
% spins_writer
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 12-Oct-2021; Last revision: 12-Oct-2021
% MATLAB Version: 9.10.0.1739362 (R2021a) Update 5

params = spins_params; 
if nargin<3
    save = false;
end

switch lower(derivative)
    case 'ri' % Richardson Number
        g_rho0 = -params.g/params.rho_0;
        u = spins_reader_new('u', ii);        
        rho = (spins_reader_new('rho', ii));
        [~, du_dz] = get_grad2(u);
        [~, drho_dz] = get_grad2(rho);
        N_sq = g_rho0 * drho_dz;
        data = N_sq./(du_dz.^2);
        
    case 'vorty' % Vorticity
        u = spins_reader_new('u', ii);        
        w = spins_reader_new('w', ii);
        [~, du_dz] = get_grad2(u);
        [dw_dx] = get_grad2(w);
        data = du_dz-dw_dx;
        
    case 'diss' % Dissipation
        nu = params.visco*params.rho_0;
        u = spins_reader_new('u', ii);
        w = spins_reader_new('w', ii);
        [du_dx, du_dz] = get_grad2(u);
        [dw_dx, dw_dz] = get_grad2(w);
        data = 2*nu*(du_dx.^2 + dw_dz.^2 + 2*(.5*(du_dz + dw_dx)).^2);
    otherwise
        split_string = split(lower(derivative), '_');
        if ~(length(split_string)==2)
            error('derivative type not configured, try ri, vorty, diss, or var_x type')
        end

        raw_data = spins_reader_new(split_string{1}, ii);
        
        [~, z] = spinsgrid2d;
%        raw_data = raw_data + -0.5*params.delta_rho * tanh((z-params.rho_loc)/params.dz_rho);
        switch split_string{2}
            case 'x'
                data = get_grad2(raw_data);
            case 'z'
                [~, data] = get_grad2(raw_data);
            case 'zz'
                [~, data] = get_grad2(raw_data);
                [~, data] = get_grad2(data);
            case 'absgrad'
                [d_dx, d_dz] = get_grad2(raw_data);
                data = sqrt(d_dx.^2 + d_dz.^2);
            otherwise
                error('derivative type not configured, try ri, vorty, diss, or var_x type')

        end
end
% Save as SPINS output type file
filename = [derivative,'.', num2str(ii)];
if save
    spins_writer(filename, data);
end