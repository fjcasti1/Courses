function getMovieFrame(u,v,Y,n,path,time)
axisSize=14;
markersize=6;
linewidth=3.5;

% Define the points of the different meshes
global xu;
global yu;
global xv;
global yv;
global xY;
global yY;

% Plot u
h1 = figure();set(gcf,'Visible', 'off');
pcolor(xu,yu,u')
shading interp
hold on
% Inlet left
plot(0:0.001:0.06,0.5,'k.','markersize',markersize,'linewidth',linewidth)
plot(0:0.001:0.06,1,'k.','markersize',markersize,'linewidth',linewidth)
% Inlet right
plot(3.94:0.001:4,1,'k.','markersize',markersize,'linewidth',linewidth)
plot(3.94:0.001:4,1.5,'k.','markersize',markersize,'linewidth',linewidth)
% Inlet top
plot(0.5,1.94:0.001:2,'k.','markersize',markersize,'linewidth',linewidth)
plot(1,1.94:0.001:2,'k.','markersize',markersize,'linewidth',linewidth)
% outlet
plot(1.5,0:0.001:0.06,'k.','markersize',markersize,'linewidth',linewidth)
plot(2.5,0:0.001:0.06,'k.','markersize',markersize,'linewidth',linewidth)
axis([0 4 0 2])
caxis([-1 1])
colorbar
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title(['t=',num2str(time)])
% yticks([0 0.25 0.50 0.75 1.0 1.25 1.50 1.75 2.0]);
% xticks([0 0.50 1.0 1.50 2.0 2.50 3.0 3.50 4.0]);
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',axisSize)
pbaspect([2 1 1])
txt=[path,'u_%03d.png'];
saveas(gcf,sprintf(txt,n))

% Plot v
h2 = figure();set(gcf,'Visible', 'off');
pcolor(xv,yv,v')
shading interp
hold on
% Inlet left
plot(0:0.001:0.06,0.5,'k.','markersize',markersize,'linewidth',linewidth)
plot(0:0.001:0.06,1,'k.','markersize',markersize,'linewidth',linewidth)
% Inlet right
plot(3.94:0.001:4,1,'k.','markersize',markersize,'linewidth',linewidth)
plot(3.94:0.001:4,1.5,'k.','markersize',markersize,'linewidth',linewidth)
% Inlet top
plot(0.5,1.94:0.001:2,'k.','markersize',markersize,'linewidth',linewidth)
plot(1,1.94:0.001:2,'k.','markersize',markersize,'linewidth',linewidth)
% outlet
plot(1.5,0:0.001:0.06,'k.','markersize',markersize,'linewidth',linewidth)
plot(2.5,0:0.001:0.06,'k.','markersize',markersize,'linewidth',linewidth)
axis([0 4 0 2])
caxis([-1.5 0.2])
colorbar
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title(['t=',num2str(time)])
% yticks([0 0.25 0.50 0.75 1.0 1.25 1.50 1.75 2.0]);
% xticks([0 0.50 1.0 1.50 2.0 2.50 3.0 3.50 4.0]);
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',axisSize)
pbaspect([2 1 1])
txt=[path,'v_%03d.png'];
saveas(gcf,sprintf(txt,n))

% Plot Y
h3 = figure();set(gcf,'Visible', 'off');
pcolor(xY,yY,Y')
shading interp
hold on
% Inlet left
plot(0:0.001:0.06,0.5,'k.','markersize',markersize,'linewidth',linewidth)
plot(0:0.001:0.06,1,'k.','markersize',markersize,'linewidth',linewidth)
% Inlet right
plot(3.94:0.001:4,1,'k.','markersize',markersize,'linewidth',linewidth)
plot(3.94:0.001:4,1.5,'k.','markersize',markersize,'linewidth',linewidth)
% Inlet top
plot(0.5,1.94:0.001:2,'k.','markersize',markersize,'linewidth',linewidth)
plot(1,1.94:0.001:2,'k.','markersize',markersize,'linewidth',linewidth)
% outlet
plot(1.5,0:0.001:0.06,'k.','markersize',markersize,'linewidth',linewidth)
plot(2.5,0:0.001:0.06,'k.','markersize',markersize,'linewidth',linewidth)
axis([0 4 0 2])
caxis([0 1])
colorbar
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title(['t=',num2str(time)])
% yticks([0 0.25 0.50 0.75 1.0 1.25 1.50 1.75 2.0]);
% xticks([0 0.50 1.0 1.50 2.0 2.50 3.0 3.50 4.0]);
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',axisSize)
pbaspect([2 1 1])
txt=[path,'Y_%03d.png'];
saveas(gcf,sprintf(txt,n))

% Plot Ymix
h4 = figure();set(gcf,'Visible', 'off');
mix=Y.*(1-Y);
pcolor(xY,yY,mix')
shading interp
hold on
% Inlet left
plot(0:0.001:0.06,0.5,'k.','markersize',markersize,'linewidth',linewidth)
plot(0:0.001:0.06,1,'k.','markersize',markersize,'linewidth',linewidth)
% Inlet right
plot(3.94:0.001:4,1,'k.','markersize',markersize,'linewidth',linewidth)
plot(3.94:0.001:4,1.5,'k.','markersize',markersize,'linewidth',linewidth)
% Inlet top
plot(0.5,1.94:0.001:2,'k.','markersize',markersize,'linewidth',linewidth)
plot(1,1.94:0.001:2,'k.','markersize',markersize,'linewidth',linewidth)
% outlet
plot(1.5,0:0.001:0.06,'k.','markersize',markersize,'linewidth',linewidth)
plot(2.5,0:0.001:0.06,'k.','markersize',markersize,'linewidth',linewidth)
axis([0 4 0 2])
caxis([0 0.25])
colorbar
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title(['t=',num2str(time)])
% yticks([0 0.25 0.50 0.75 1.0 1.25 1.50 1.75 2.0]);
% xticks([0 0.50 1.0 1.50 2.0 2.50 3.0 3.50 4.0]);
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',axisSize)
pbaspect([2 1 1])
txt=[path,'Ymix_%03d.png'];
saveas(gcf,sprintf(txt,n))
end