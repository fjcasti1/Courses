function ContourPlots(u,v,Y,M,N,hx,hy,n)
axisSize=14;
markersize=16;
linewidth=3.5;

% Define the points of the different meshes
xu=linspace(0,4,M+1);
yu=linspace(-hy/2,2+hy/2,N+2);
xv=linspace(-hx/2,4+hx/2,M+2);
yv=linspace(0,2,N+1);
xY=linspace(-5*hx/2,4+5*hx/2,M+6);
yY=linspace(-5*hx/2,2+5*hx/2,N+6);

% Plot u
figure(10+n)
contourf(xu,yu,u')
hold on
plot(1,0.5,'r.','markersize',markersize,'linewidth',linewidth)
axis([0 4 0 2])
caxis([-1.5 3])
colorbar
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
% yticks([0 0.25 0.50 0.75 1.0 1.25 1.50 1.75 2.0]);
% xticks([0 0.50 1.0 1.50 2.0 2.50 3.0 3.50 4.0]);
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',axisSize)
pbaspect([2 1 1])
grid on
txt=['Latex/FIGURES/u_' num2str(10+n)];
saveas(gcf,txt,'epsc')

% Plot v
figure(20+n)
contourf(xv,yv,v')
hold on
plot(1,1.5,'r.','markersize',markersize,'linewidth',linewidth)
axis([0 4 0 2])
caxis([-1.5 0])
colorbar
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
% yticks([0 0.25 0.50 0.75 1.0 1.25 1.50 1.75 2.0]);
% xticks([0 0.50 1.0 1.50 2.0 2.50 3.0 3.50 4.0]);
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',axisSize)
pbaspect([2 1 1])
grid on
txt=['Latex/FIGURES/v_' num2str(20+n)];
saveas(gcf,txt,'epsc')

% Plot Y
figure(30+n)
contourf(xu,yu,u')
hold on
plot(2,0.5,'r.','markersize',markersize,'linewidth',linewidth)
axis([0 4 0 2])
caxis([-0.1 1.1])
colorbar
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
% yticks([0 0.25 0.50 0.75 1.0 1.25 1.50 1.75 2.0]);
% xticks([0 0.50 1.0 1.50 2.0 2.50 3.0 3.50 4.0]);
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',axisSize)
pbaspect([2 1 1])
grid on
txt=['Latex/FIGURES/Y_' num2str(30+n)];
saveas(gcf,txt,'epsc')
end