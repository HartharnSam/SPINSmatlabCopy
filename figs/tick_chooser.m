function tick_chooser(axesticks, ax)
%TICK_CHOOSER - Makes nice + logically spaced ticks on the chosen axis
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
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 14-Dec-2020; Last revision: 14-Dec-2020
% MATLAB Version: 9.9.0.1467703 (R2020b)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
if nargin<1
    ax = gca;
end
if strcmpi(axesticks, 'XTick')
    limits = ax.XLim;
    width = ax.Position(3);
elseif strcmpi(axesticks, 'YTick')
    limits = ax.YLim;
    width = ax.Position(3);
else
    error("Axesticks must be 'XTick' or 'YTick'")
end

% Identify ideal number of ticks
if width < .3
    low_limit = 3;
else
    low_limit = 6; 
end

upp_limit = low_limit*1.333; 
% First approximation of tick spacing
tick_space = (diff(limits)/low_limit); % Assumes there can be 6 ticks
% Choose a tick spacing from the appropriate list
reasonable_diffs = [0.1 .2 .5 1 2 5 10 15 20 25 50 100 150 200 250 500];
[~, tick_min_ind] = sort(abs(tick_space - reasonable_diffs));
tick_space = reasonable_diffs(tick_min_ind(1));
% Check tick spacing works, iterate if not
i = 2;
while diff(limits)/tick_space > upp_limit
    tick_space = reasonable_diffs(tick_min_ind(i));
    i = i+1;
end

% Set tick locations
tick_locs = round((limits(1):tick_space:limits(2))/tick_space)*tick_space;
set(ax, axesticks, tick_locs);

end
