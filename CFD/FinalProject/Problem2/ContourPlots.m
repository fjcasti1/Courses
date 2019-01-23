function ContourPlots(u,v,Y,n,path)
axisSize=14;
markersize=6;
linewidth=3.5;

w=0.06;

% Define the points of the different meshes
global Lx; global Ly;
global xu; global yu;
global xv; global yv;
global xY; global yY;
global a1; global b1;
global a2; global b2;
global a3; global b3;
global ao; global bo;

% Plot u
figure();set(gcf,'Visible', 'off');
pcolor(xu,yu,u')
shading interp
hold on
% Inlet left
plot(0:0.001:w,a1,'k.','markersize',markersize,'linewidth',linewidth)
plot(0:0.001:w,b1,'k.','markersize',markersize,'linewidth',linewidth)
% Inlet right
plot(Lx-w:0.001:Lx,a2,'k.','markersize',markersize,'linewidth',linewidth)
plot(Lx-w:0.001:Lx,b2,'k.','markersize',markersize,'linewidth',linewidth)
% Inlet top
plot(a3,Ly-w:0.001:Ly,'k.','markersize',markersize,'linewidth',linewidth)
plot(b3,Ly-w:0.001:Ly,'k.','markersize',markersize,'linewidth',linewidth)
% outlet
plot(ao,0:0.001:w,'k.','markersize',markersize,'linewidth',linewidth)
plot(bo,0:0.001:w,'k.','markersize',markersize,'linewidth',linewidth)
axis([0 Lx 0 Ly])
caxis([-3 3])
colorbar
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
% yticks([0 0.25 0.50 0.75 1.0 1.25 1.50 1.75 2.0]);
% xticks([0 0.50 1.0 1.50 2.0 2.50 3.0 3.50 4.0]);
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',axisSize)
pbaspect([1 2 1])
% grid on
txt=[path,'u_' num2str(n)];
saveas(gcf,txt,'epsc')

% Plot v
figure();set(gcf,'Visible', 'off');
pcolor(xv,yv,v')
shading interp
hold on
% Inlet left
plot(0:0.001:w,a1,'k.','markersize',markersize,'linewidth',linewidth)
plot(0:0.001:w,b1,'k.','markersize',markersize,'linewidth',linewidth)
% Inlet right
plot(Lx-w:0.001:Lx,a2,'k.','markersize',markersize,'linewidth',linewidth)
plot(Lx-w:0.001:Lx,b2,'k.','markersize',markersize,'linewidth',linewidth)
% Inlet top
plot(a3,Ly-w:0.001:Ly,'k.','markersize',markersize,'linewidth',linewidth)
plot(b3,Ly-w:0.001:Ly,'k.','markersize',markersize,'linewidth',linewidth)
% outlet
plot(ao,0:0.001:w,'k.','markersize',markersize,'linewidth',linewidth)
plot(bo,0:0.001:w,'k.','markersize',markersize,'linewidth',linewidth)
axis([0 Lx 0 Ly])
caxis([-3 0.35])
colorbar
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
% yticks([0 0.25 0.50 0.75 1.0 1.25 1.50 1.75 2.0]);
% xticks([0 0.50 1.0 1.50 2.0 2.50 3.0 3.50 4.0]);
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',axisSize)
pbaspect([1 2 1])
% grid on
txt=[path,'v_' num2str(n)];
saveas(gcf,txt,'epsc')

% Plot Y
figure();set(gcf,'Visible', 'off');
pcolor(xY,yY,Y')
shading interp
hold on
% Inlet left
plot(0:0.001:w,a1,'k.','markersize',markersize,'linewidth',linewidth)
plot(0:0.001:w,b1,'k.','markersize',markersize,'linewidth',linewidth)
% Inlet right
plot(Lx-w:0.001:Lx,a2,'k.','markersize',markersize,'linewidth',linewidth)
plot(Lx-w:0.001:Lx,b2,'k.','markersize',markersize,'linewidth',linewidth)
% Inlet top
plot(a3,Ly-w:0.001:Ly,'k.','markersize',markersize,'linewidth',linewidth)
plot(b3,Ly-w:0.001:Ly,'k.','markersize',markersize,'linewidth',linewidth)
% outlet
plot(ao,0:0.001:w,'k.','markersize',markersize,'linewidth',linewidth)
plot(bo,0:0.001:w,'k.','markersize',markersize,'linewidth',linewidth)
axis([0 Lx 0 Ly])
caxis([0 1])
colorbar
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
% yticks([0 0.25 0.50 0.75 1.0 1.25 1.50 1.75 2.0]);
% xticks([0 0.50 1.0 1.50 2.0 2.50 3.0 3.50 4.0]);
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',axisSize)
pbaspect([1 2 1])
% grid on
txt=[path,'Y_' num2str(n)];
saveas(gcf,txt,'epsc')

% Plot Ymix
figure();set(gcf,'Visible', 'off');
mix=Y.*(1-Y);
pcolor(xY,yY,mix')
shading interp
hold on
% Inlet left
plot(0:0.001:w,a1,'k.','markersize',markersize,'linewidth',linewidth)
plot(0:0.001:w,b1,'k.','markersize',markersize,'linewidth',linewidth)
% Inlet right
plot(Lx-w:0.001:Lx,a2,'k.','markersize',markersize,'linewidth',linewidth)
plot(Lx-w:0.001:Lx,b2,'k.','markersize',markersize,'linewidth',linewidth)
% Inlet top
plot(a3,Ly-w:0.001:Ly,'k.','markersize',markersize,'linewidth',linewidth)
plot(b3,Ly-w:0.001:Ly,'k.','markersize',markersize,'linewidth',linewidth)
% outlet
plot(ao,0:0.001:w,'k.','markersize',markersize,'linewidth',linewidth)
plot(bo,0:0.001:w,'k.','markersize',markersize,'linewidth',linewidth)
axis([0 Lx 0 Ly])
caxis([0 0.25])
colorbar
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
% yticks([0 0.25 0.50 0.75 1.0 1.25 1.50 1.75 2.0]);
% xticks([0 0.50 1.0 1.50 2.0 2.50 3.0 3.50 4.0]);
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',axisSize)
pbaspect([1 2 1])
txt=[path,'Ymix_' num2str(n)];
saveas(gcf,txt,'epsc')
end