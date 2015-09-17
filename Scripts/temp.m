clc
clear 
close all

% 
% lvl_figure_handle = figure('Position', window_pos4, ...
%         'Renderer', 'OpenGL', 'MenuBar', 'none', 'ToolBar', 'none');
    
    % Add frequency spectrum plot to the figure
%     lvl_axes_handle = axes('Parent', lvl_figure_handle);
    m=[5 0; 0 5];
  bar_handle = createfigure(m);
   ylim([1 10]) 
    title('attention and meditation level','FontSize', 18);
    
    ylabel('level', 'FontSize', 16);
    i=0;
    while 1
        
        i=i+1;
        for k=1:length(bar_handle)
        mm=[rand(1)*10 0; 0 rand(1)*10];
        
        set(bar_handle(k),'Ydata', mm(k,:));
%         get(bar_handle(k),'Ydata')
        end
        drawnow expose
    end