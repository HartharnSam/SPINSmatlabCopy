function data = SPINS_derivs(derivative, ii, save)
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

switch lower(derivative)
    case 'ri' % Richardson Number
        g_rho0 = -params.g/params.rho_0;
                u = spins_reader_new('u', ii);        
                rho = spins_reader_new('rho', ii);
        [~, du_dz] = get_grad2(u);
        [~, drho_dz] = get_grad2(rho);
        data = g_rho0*drho_dz.*(du_dz.^2);
        
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

end
% Save as SPINS output type file
filename = [derivative,'.', num2str(ii)];
if save
    spins_writer(filename, data);
end