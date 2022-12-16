function [ax1, clrbar_ax] = plot_ri(x, z, data)

%Testdata
testdata = -.5:0.05:2;
testdata = repmat(testdata, 10, 1);
% data = testdata;
% x = testdata(1, :)
% z = 1:10;
levels = [0 .25 1];

data(data<levels(1)) = -.5; 
testdata(testdata<levels(1)) = -0.5;

data(data>levels(end)) = 1.5;
testdata(testdata>levels(end)) = 1.5;

for ii = 1:length(levels)-1
    data(data >= levels(ii) & data < levels(ii+1)) = (levels(ii)+levels(ii+1))/2;
    testdata(testdata >= levels(ii) & testdata < levels(ii+1)) = (levels(ii)+levels(ii+1))/2;

end

%pcolor(-.5:0.01:2, [1:10], data);
[cf cfh] = contourf(x, z, data, levels);
cfh.LineStyle = 'none';
ax1= gca;
clrmap = custom_color_orders('set_colormap', 'Art', 'Darjeeling1',5)
clrmap = clrmap(1:4, :);
colormap(gca, clrmap);
c = colorbar;
c_loc = c.Position;
ax_pos = ax1.Position;

delete(c);

clrbar_ax = axes;
ax1.Position = ax_pos;
clrbar_ax.Position = c_loc;
clrmap = colormap(gca, clrmap);

contourf([1:10], -.5:0.05:2, testdata', levels);
xticks([]); clrbar_ax.YAxisLocation = 'right';
yticks([-.25 levels(1:end-1)+diff(levels)/2 1.5]);
yticklabels({'<0', '0-0.25', '0.25 - 1', '>1'});
ylabel('$Ri$', 'interpreter', 'latex');

axes(ax1);