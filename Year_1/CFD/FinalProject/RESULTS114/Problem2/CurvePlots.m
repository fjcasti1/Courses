function CurvePlots(R,S,t,path)
axisSize=14;
linewidth=2;

% Cut out what hasn't been filled
index=find(t==0,2);
index=index(2);
t(index:end)=[];
R(index:end)=[];
S(index:end)=[];

% Performance parameter R
figure();set(gcf,'Visible', 'off');
plot(t,R,'linewidth',linewidth)
xlim([0 t(end)])
xlabel('$t$','Interpreter','latex')
ylabel('$R(t)$','Interpreter','latex')
set(gca,'fontsize',axisSize)
grid on
txt=[path,'Rplot'];
saveas(gcf,txt,'epsc')

% Performance parameter S
figure();set(gcf,'Visible', 'off');
plot(t,S,'linewidth',linewidth)
xlim([0 t(end)])
xlabel('$t$','Interpreter','latex')
ylabel('$S(t)$','Interpreter','latex')
set(gca,'fontsize',axisSize)
grid on
txt=[path,'Splot'];
saveas(gcf,txt,'epsc')
end