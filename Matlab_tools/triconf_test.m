%% Read in data
%[x,y,d,mesh] = read_in_pde_vtk('filename','MY_psi_back_info.vtk');
[x,y,d,mesh] = read_in_pde_vtk('filename','vtkFIG3demo.vtk');

%% Manipulate data
close all
figure(1)
axis equal
hold on
%[CS,h]=tricontf(x,y,mesh,d,8-[2.254,2.71129,2.85493,3.00003,5.746]);
%[CS,h]=tricontf(x,y,mesh,d,[max(d)/4,max(d)/2,3*max(d)/4]);
%[CS,h]=tricontf(x,y,mesh,d);
[CS,h]=tricontf(x,y,mesh,d,linspace(1,99,8)/100*max(d));
set(h,'edgecolor','none')
set(gcf,'Renderer','painters')

% trimesh(mesh,x,y,'Color','k','LineWidth',1);

% To get individual contour levels as lines, get XData and YData from the
% appropriate h patch. This includes boundary info. Can delete this by
% comparing to patch 1 like so:

p1 = patch('Faces',mesh,'Vertices',[x,y]);
set(p1,'FaceColor','none','EdgeColor','black')

%% Debugging
figure(2)
axis equal
hold on

 h_1 = h(1);
 h_2 = h(2);
for j = 1:numel(h)

% 
% BoundaryData = [h_num.XData,h_1.YData];
% ContourData = [h_num.XData,h_num.YData];
% 
% FinalContour = setdiff(ContourData,BoundaryData,'rows');
% plot(FinalContour(:,1),FinalContour(:,2),'r');

plot(h(j).XData,h(j).YData,'r');

pause

end