# SPINSmatlab
Project for plotting and analyzing SPINS outputs in matlab

Usage:
First read in grid:
  gdpar = spins_gridparams();         % for unmapped grid
  gdpar = spins_gridparams('Full');   % for mapped grid
(No arguments will give vectorized grid. This is useful for unmapped grids. Mapped grids require the full grid to be loaded)

Second make plot:
  spins_plot2d(var, t_i);   % plot var at output t_i, see file for more documentation
