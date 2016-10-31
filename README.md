# SPINSmatlab
Project for plotting and analyzing SPINS outputs in matlab

Usage:

    1. Read in grid:

| Command                               | Grid |
| ---                                   | --- |
| gdpar = spins_gridparams();           | unmapped |  
| gdpar = spins_gridparams('Vector');   | unmapped (default) |  
| gdpar = spins_gridparams('Full');		  | mapped |  
| gdpar = spins_gridparams('FastFull');	| mapped |  

Passing no argument will default to type 'Vector'. This is useful for unmapped grids which only require 1D vectors.  
Mapped grids require the full grid to be loaded. 'Full' reads the entire grids.  
2D plotting only requires cross sections, so 'FastFull' searches the spins.conf for parameters for the x-y slices. The x-z and y-z slices are read from the grids. Unless you require the full field, 'FastFull' will be sufficient.

    2. Make plot:  
| Command                     | Purpose |
| ---                         | --- |
| spins_plot2d(var, t_i);     | plot var at output number t_i |  
| spins_plot2d(var, 't');     | plot var at time t |  

There are many optional arguments for adjusting the plot, see spins_plot2d.m for more documentation.
