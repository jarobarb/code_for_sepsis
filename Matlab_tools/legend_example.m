close(figure(1)); figure(1);
x = linspace(0,2*pi);
y1 = sin(x);
y2 = cos(x);
y3 = (x-pi)./2;
plot(x,y1,x,y2,x,y3);
pixelsperinch = 128;
set(gcf,'Position',pixelsperinch*[1,1,6,4]);

% legend('y1','y2','y3');
% %  Gets "handle" or numeric label associated with the axis
% axis_handle = gca;
% 
% %  Gets "handle" or numeric label associated with the lines
% line_handles = get(gca,'Children');
% 
% %  Turn off the legend for y2
% line_handles(2).Annotation.LegendInformation.IconDisplayStyle = 'off';