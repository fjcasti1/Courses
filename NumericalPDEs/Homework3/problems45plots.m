clear all
close all
clc
linewidth=1.6;
labelfontsize=18;
legendfontsize=14;

sol0=load('dat/sol0.dat');
sol05=load('dat/sol05.dat');
sol1=load('dat/sol1.dat');
sol15=load('dat/sol15.dat');
sol2=load('dat/sol2.dat');

solDx25=load('dat/solDx25.dat');
solDx50=load('dat/solDx50.dat');
solDx100=load('dat/solDx100.dat');
solDx200=load('dat/solDx200.dat');

figure(1)
plot(sol0(:,1),sol0(:,2),'linewidth',linewidth)
hold on
plot(sol05(:,1),sol05(:,2),'linewidth',linewidth)
plot(sol1(:,1),sol1(:,2),'linewidth',linewidth)
plot(sol15(:,1),sol15(:,2),'linewidth',linewidth)
plot(sol2(:,1),sol2(:,2),'linewidth',linewidth)
xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
ylabel('$u(x,t)$','interpreter','latex','fontsize',labelfontsize)
grid on
h1=legend('$t_{comp}=0$','$t_{comp}=0.5$','$t_{comp}=1$','$t_{comp}=1.5$','$t_{comp}=2$');
set(h1,'interpreter','latex','fontsize',legendfontsize)
saveas(gcf,'problem4','png')

figure(2)
plot(solDx25(:,1),solDx25(:,2),'linewidth',linewidth)
hold on
plot(solDx50(:,1),solDx50(:,2),'linewidth',linewidth)
plot(solDx100(:,1),solDx100(:,2),'linewidth',linewidth)
plot(solDx200(:,1),solDx200(:,2),'linewidth',linewidth)
xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
ylabel('$u(x,t)$','interpreter','latex','fontsize',labelfontsize)
grid on
h1=legend('$N=25$','$N=50$','$N=100$','$N=200$');
set(h1,'interpreter','latex','fontsize',legendfontsize)
saveas(gcf,'problem5','png')

figure(3)
plot(solDx25(:,1),solDx25(:,2),'linewidth',linewidth)
hold on
plot(solDx50(:,1),solDx50(:,2),'linewidth',linewidth)
plot(solDx100(:,1),solDx100(:,2),'linewidth',linewidth)
plot(solDx200(:,1),solDx200(:,2),'linewidth',linewidth)
xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
ylabel('$u(x,t)$','interpreter','latex','fontsize',labelfontsize)
grid on
axis([0 4 0 0.9])
h1=legend('$N=25$','$N=50$','$N=100$','$N=200$');
set(h1,'interpreter','latex','fontsize',legendfontsize)
saveas(gcf,'problem5Zoom','png')