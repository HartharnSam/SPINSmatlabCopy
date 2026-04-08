function plot_khbillow(ii, ylims, doRi)
%% PLOT_KHBILLOW_CASE
% Plotting tool for the kh_billow case where rho output is a rho', and
% needs to be re-produced
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% Department of Geography & Environmental Sciences, Northumbria University
% email address: sam.hartharn-evans@northumbria.ac.uk
% GitHub: https://github.com/HartharnSam

clc; close all;
spinsstartup;

[x, z] = spinsgrid2d;
times = get_output_times;
params = spins_params;

vorty = spins_derivs('vorty', ii);
rho = spins_reader_new('rho', ii);
rho = rho + -0.5*params.delta_rho * tanh((z-params.rho_loc)/params.dz_rho);

%rho_z = spins_derivs('rho_z', ii);

u = spins_reader_new('u', ii);

if doRi
    ri = spins_derivs('Ri', ii);
    ri_mean = mean(ri);
    ri_mean_nonan_i = find(~isnan(ri_mean) & ri_mean>0);
    ri_mean_nonan = ri_mean(ri_mean_nonan_i);

    zci = @(v) find(diff(sign(v)));
    regions = unique([1,  zci(ri_mean_nonan-.25), length(ri_mean_nonan)]);
    region_sign = sign(ri_mean_nonan(regions(1:end-1)+1)-.25);
    phys_regions = z(1, ri_mean_nonan_i(regions));
end

figure(1);
subplot(1, 3, 1)
plot(u(1, :), z(1, :));
hold on;
title("t = "+times(ii+1))

if doRi
    for j = 1:length(phys_regions)-1
        if region_sign(j)<0
            p = fill([min(xlim) max(xlim) max(xlim) min(xlim)], [phys_regions(j) phys_regions(j) phys_regions(j+1) phys_regions(j+1)], [0.8 0.8 0.8]);
            p.FaceAlpha = .4; p.EdgeColor = 'none';
        end
    end
end
xlim([min(xlim) max(xlim)])
ylim(ylims)

subplot(1, 3, 2); hold on;
plot(rho(1, :), z(1, :));
ylim(ylims)

if doRi
subplot(1, 3, 3); hold on;
plot(ri_mean, z(1, :))
xline(.25)
xlim([0 1])
ylim(ylims)

    for j = 1:length(phys_regions)-1
        if region_sign(j)<0
            p = fill([min(xlim) max(xlim) max(xlim) min(xlim)], [phys_regions(j) phys_regions(j) phys_regions(j+1) phys_regions(j+1)], [0.8 0.8 0.8]);
            p.FaceAlpha = .4; p.EdgeColor = 'none';
        end
    end
end
figure(2)
subplot(2,1, 1)
pcolor(x, z, rho);
ylim(ylims)
params = spins_params;
%rho_piv = (params.rho_top+params.rho_bot)/2;
cmocean('balance');

title("t = "+times(ii+1))

subplot(2,1,2);
%w = spins_reader_new('w', ii);
pcolor(x, z, vorty);
ylim(ylims)
%clim([-5 25])
cmocean('balance', 'pivot', 0);
%colorbar; 