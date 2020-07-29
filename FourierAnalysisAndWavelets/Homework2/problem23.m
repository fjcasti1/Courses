linewidth=2;
fontsize=15;
x=linspace(-3*pi,3*pi,1000);
y=cos(x)+abs(cos(x));
plot(x,y,'linewidth',linewidth)
axis([-3*pi 3*pi -1 3])
grid on
xlabel('$x$','interpreter','latex','fontsize',fontsize)
ylabel('$f(x)$','interpreter','latex','fontsize',fontsize)
h=legend('$f(x)=\cos(x)+|\cos(x)|$');
set(h,'interpreter','latex')
saveas(gcf,'problem23','epsc')