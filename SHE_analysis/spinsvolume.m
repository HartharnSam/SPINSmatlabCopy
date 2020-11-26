function volume = spinsvolume()
%FUNCTION_NAME - One line description of what the function or script performs (H1 line)%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    input1 - Description
%    input2 - Description
%    input3 - Description
%
% Outputs:
%    output1 - Description
%    output2 - Description
%
% Example:
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Denis Gilbert, Ph.D., physical oceanography
% Maurice Lamontagne Institute, Dept. of Fisheries and Oceans Canada
% email address: gilbertd@dfo-mpo.gc.ca
% GitHub: http://www.qc.dfo-mpo.gc.ca/iml/
% December 1999; Last revision: 12-May-2004
% MATLAB Version: Release: Service Pack

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
params = spins_params;
spinsgrid2d;
xi = [zeros(1, NZ); x; ones(1, NZ)*Lx];
dx = diff(xi);
zi = [ones(NX, 1)*-Lz z zeros(NX, 1)];
dz = diff(zi, 1, 2);

dx2 = dx/2;
dz2 = dz/2;

volume = (dx2(1:end-1, :)+dx2(2:end, :)).*(dz2(:, 1:end-1)+dz2(:, 2:end));

%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------