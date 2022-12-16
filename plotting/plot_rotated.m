% Plot rotated by slope angle
clc; clearvars; close all; 

spins_gridparams('full');

global gdpar;
gdpar.params.tilt_angle = atand(gdpar.params.hill_slope);
gdpar_rot = rotate_refframe;

rho = spins_reader_new('rho', 150);
figure
subplot(2, 1, 1)
pcolor(gdpar_rot.gd.x, gdpar_rot.gd.z, rho); shading flat;
xlim([5.4 14.4])
ylim([0 .3])

subplot(2, 1, 2)
pcolor(gdpar.gd.x, gdpar.gd.z, rho);
shading flat;
xlim([5.5 14.5])
ylim([-0.3 0])