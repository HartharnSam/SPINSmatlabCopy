function settings = betterplots(target)
settings = struct('DefaultLineLineWidth',2,'DefaultTextFontSize',12,...
    'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
    'DefaultAxesFontWeight','bold', 'defaultfigurecolor', 'w',...
    'DefaultAxesXDir', 'reverse', 'DefaultFigureColormap', darkjet,...
    'DefaultAxesBox', 'on', 'DefaultAxesNextPlot', 'replacechildren');

if nargin<1
    target = gcf;
end
fields = fieldnames(settings);
for i = 1:length(fields)
    str = join(['settings.', fields(i),';'], '');
    set(target, fields{i}, eval(str{1}));
end
