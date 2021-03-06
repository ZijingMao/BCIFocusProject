function [bar1] = createfigure(ymatrix1)
%CREATEFIGURE(YMATRIX1)
%  YMATRIX1:  bar matrix data

%  Auto-generated by MATLAB on 24-Jul-2015 16:58:05

% Create figure
% figure1 = figure;
window_pos4 = [800, 100, 500, 400];
figure1 = figure('Position', window_pos4, ...
        'Renderer', 'OpenGL', 'MenuBar', 'none', 'ToolBar', 'none');
% Create axes
axes1 = axes('Parent',figure1,'XTick',[1 2]);
%% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 3]);
%% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0 5]);
box(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to bar
bar1 = bar(ymatrix1,'BarLayout','stacked','Parent',axes1);
set(bar1(1),'DisplayName','Attention');
set(bar1(2),'DisplayName','Meditation');
ylim([0, 5]);
xlim([0, 3]);

% Create legend
legend(axes1,'show');

