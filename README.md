# SPINSmatlab
Project for plotting and analyzing SPINS outputs in matlab

Usage:

    1. Read in grid:

| Command                               | Grid |
| ---                                   | --- |
| gdpar = spins_gridparams();           | unmapped |  
| gdpar = spins_gridparams('Vector');   | unmapped (default) |  
| gdpar = spins_gridparams('Full');		| mapped |  
| gdpar = spins_gridparams('FastFull');	| mapped |  

Passing no argument will default to type 'Vector'. This is useful for unmapped grids which only require 1D vectors.  
Mapped grids require the full grid to be loaded. 'Full' reads the entire grids.  
'FastFull' also creates a full grid by creating the grid from the parameters in spins.conf rather than reading the grids. This works for the x and y grids which are topography invariant. The z grid is still created by reading an x-z slice and using repmat to fill out the 3rd dimension.

    2. Make plot:  
| Command                     | Purpose |
| ---                         | --- |
| spins_plot2d(var, t_i);     | plot var at output number t_i |  
| spins_plot2d(var, 't');     | plot var at time t |  

There are many optional arguments for adjusting the plot, see spins_plot2d.m for more documentation.

To automatically include the SPINSmatlab directory in the MATLAB search path, add the following to `startup.m` which must exist in the MATLAB path.

```
% add path to SPINSmatlab
sm_loc = genpath('/path/to/SPINSmatlab');
addpath(sm_loc);
% remove all hidden subdirectories
all_subs =  strsplit(sm_loc,':');
for ii = 1:length(all_subs)
    if ~isempty(strfind(all_subs{ii},'.'))
        rmpath(all_subs{ii});
    end
end
% clean work space
clear all
```
