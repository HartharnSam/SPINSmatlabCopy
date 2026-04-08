function spins_clean
%SPINS_CLEAN - Automatic cleanup function for folders - removing actual
%rundata and keeping only the spins.conf, .cpp etc type files
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% Department of Geography & Environmental Sciences, Northumbria University
% email address: sam.hartharn-evans@northumbria.ac.uk
% GitHub: https://github.com/HartharnSam
% 05-Aug-2024; Last revision: 05-Aug-2024
% MATLAB Version: 9.10.0.1739362 (R2021a) Update 5

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
delete w.* u.* v.* rho.* t.* s.* % delete variable outputs
delete vorty.* diss.* % delete derivatives outputs
delete *.txt *.out % delete text based outputs
delete *grid* % delete grid outputs
%delete *.mp4 *.png % delete processed data
delete *.mat % and more processed data
delete *.x *.src.* *.o % delete any executables
%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------