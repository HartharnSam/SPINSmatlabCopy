clc; clearvars; close all; 
cd C:\Users\samha\Documents\PhD\Numerics\Model\lx_test_4
ii = 38;
figure(1)
rho_0 = 1026;
delta_rho = 0.019493;

paths = {'../lx_test_6', '../lx_test_5', '../lx_test_4', '../lx_test_1', '../amp_test_1', '../lx_test_2'};
Lxs = [.3, .35, .4, .5, .6, .7];

% Lx 1
spinsgrid2d;
for i =1:length(Lxs)
cd(paths{i})
rho=spins_reader_new('rho',ii);
rho = rho_converter(rho);
u=spins_reader_new('u',ii);
umaxabs=max(abs(u(:)));
%betterplots
subaxis(2, length(Lxs), i, 1)
    pcolor(x,z,rho),shading flat
    %colormap(gca, cmocean('dense'));
    colormap darkjet
    caxis(gca, [rho_0 1046]);
    %colorbar
    set(gca, 'XDir', 'reverse');
    title(['Lx =  ', num2str(Lxs(i))])
    %contourf(x,z,rho),shading flat % This is better for printed figures
    %ylabel('z (m)')
subaxis(2, length(Lxs), i, 2)
    pcolor(x,z,u),shading flat;
    caxis([-1 1]*umaxabs)
    %colormap(gca, cmocean('delta'));
    %ylabel('z (m)')
    xlabel('x (m)')
    set(gca, 'XDir', 'reverse');
    
end
 