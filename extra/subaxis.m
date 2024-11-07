function h = subaxis(rows, cols, cellx, celly, spanx, spany, opts)
%SUBAXIS Create axes in tiled positions. (just like subplot)
%   Usage:
%      h=subaxis(rows,cols,cellno[,settings])
%      h=subaxis(rows,cols,cellx,celly[,settings])
%      h=subaxis(rows,cols,cellx,celly,spanx,spany[,settings])
%
% SETTINGS: Spacing,SpacingHoriz,SpacingVert
%           Padding,PaddingRight,PaddingLeft,PaddingTop,PaddingBottom
%           Margin,MarginRight,MarginLeft,MarginTop,MarginBottom
%           Holdaxis
%
%           all units are relative (i.e. from 0 to 1)
%
%           Abbreviations of parameters can('t currently) be used.. (Eg MR instead of MarginRight)
%           (holdaxis means that it wont delete any axes below.)
%
%
% Example:
%
%   >> subaxis(2,1,1,'SpacingVert',0,'MR',0);
%   >> imagesc(magic(3))
%   >> subaxis(2,'p',.02);
%   >> imagesc(magic(4))
%
% 2001-2014 / Aslak Grinsted  (Feel free to modify this code.)
% 2024 / Sam Hartharn-Evans (major edits to improve argument passing and
% speed)
arguments
    rows (1, 1) double
    cols(1, 1) double
    cellx
    celly = []
    spanx = []
    spany = []
    opts.Spacing = 0.05
    opts.SpacingHoriz = 0.05
    opts.SpacingVert = 0.05
    opts.Padding = 0
    opts.PaddingRight = 0
    opts.PaddingLeft = 0
    opts.PaddingTop = 0
    opts.PaddingBottom = 0
    opts.Margin = 0.1
    opts.MarginRight = 0.1
    opts.MarginLeft = 0.1
    opts.MarginTop = 0.1
    opts.MarginBottom = 0.1
    opts.Holdaxis = 0.1
end
f = gcf;

UserDataArgsOK = 0;
Args = get(f,'UserData');

if isstruct(Args)
    UserDataArgsOK = isfield(Args,'SpacingHorizontal') & isfield(Args,'Holdaxis') & isfield(Args,'rows') & isfield(Args,'cols');
end
OKToStoreArgs = isempty(Args) | UserDataArgsOK;

%Args = parseArgs(varargin,Args,{'Holdaxis'},{'Spacing' {'sh','sv'}; 'Padding' {'pl','pr','pt','pb'}; 'Margin' {'ml','mr','mt','mb'}});

if OKToStoreArgs
    set(f,'UserData',Args);
end

%Args.rows = rows;
%Args.cols = cols;

if ~isempty(cellx) && isempty(celly)
    if numel(cellx) > 1 % restore subplot(m,n,[x y]) behaviour
        [x1, y1] = ind2sub([cols rows],cellx(1)); % subplot and ind2sub count differently (column instead of row first) --> switch cols/rows
        [x2, y2] = ind2sub([cols rows],cellx(end));
    else
        x1 = mod((cellx-1),cols)+1; x2 = x1;
        y1 = floor((cellx-1)/cols)+1; y2 = y1;
    end
elseif ~isempty(cellx) && ~isempty(celly)
        x1 = cellx; x2 = x1;
        y1 = celly; y2 = y1;  
elseif ~isempty(spanx) 
    x1 = cellx; x2 = x1+spanx-1;
    y1 = celly; y2 = y1+spany-1;
else
    error('subaxis argument error')
end

cellwidth = ((1-opts.MarginLeft-opts.MarginRight)-(cols-1)*opts.SpacingHoriz)/cols;
cellheight = ((1-opts.MarginTop-opts.MarginBottom)-(rows-1)*opts.SpacingVert)/rows;
xpos1 = opts.MarginLeft+opts.PaddingLeft+cellwidth*(x1-1)+opts.SpacingHoriz*(x1-1);
xpos2 = opts.MarginLeft-opts.PaddingRight+cellwidth*x2+opts.SpacingHoriz*(x2-1);
ypos1 = opts.MarginTop+opts.PaddingTop+cellheight*(y1-1)+opts.SpacingVert*(y1-1);
ypos2 = opts.MarginTop-opts.PaddingBottom+cellheight*y2+opts.SpacingVert*(y2-1);

if opts.Holdaxis
    h = axes('position',[xpos1 1-ypos2 xpos2-xpos1 ypos2-ypos1]);
else
    h = subplot('position',[xpos1 1-ypos2 xpos2-xpos1 ypos2-ypos1]);
end

set(h,'box','on');
%h=axes('position',[x1 1-y2 x2-x1 y2-y1]);
set(h,'units',get(gcf,'defaultaxesunits'));
set(h,'tag','subaxis');
drawnow

if (nargout == 0), clear h; end

end

