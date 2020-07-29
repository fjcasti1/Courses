function CurvePlots(pu,pv,pY,R,t)
axisSize=14;
markersize=16;
linewidth=3.5;

% Cut out what hasn't been filled
index=find(t==0,2);
index=index(2);
t(index:end)=[];
pu(index:end)=[];
pv(index:end)=[];
pY(index:end)=[];
R(index:end)=[];

% Probe for u
figure(1)
plot(t,pu,'linewidth',2)
xlim([0 t(end)])
xlabel('$t$','Interpreter','latex')
ylabel('$u(1,0.5,t)$','Interpreter','latex')
set(gca,'fontsize',axisSize)
grid on
txt='Latex/FIGURES/probeu';
saveas(gcf,txt,'epsc')

% Probe for v
figure(2)
plot(t,pv,'linewidth',2)
xlim([0 t(end)])
xlabel('$t$','Interpreter','latex')
ylabel('$v(1,1.5,t)$','Interpreter','latex')
set(gca,'fontsize',axisSize)
grid on
txt='Latex/FIGURES/probev';
saveas(gcf,txt,'epsc')

% Probe for Y
figure(3)
plot(t,pY,'linewidth',2)
xlim([0 t(end)])
xlabel('$t$','Interpreter','latex')
ylabel('$Y(2,0.5,t)$','Interpreter','latex')
set(gca,'fontsize',axisSize)
grid on
txt='Latex/FIGURES/probeY';
saveas(gcf,txt,'epsc')

% Performance parameter R
figure(4)
plot(t,R,'linewidth',2)
xlim([0 t(end)])
xlabel('$t$','Interpreter','latex')
ylabel('$R(t)$','Interpreter','latex')
set(gca,'fontsize',axisSize)
grid on
txt='Latex/FIGURES/Rplot';
saveas(gcf,txt,'epsc')
end