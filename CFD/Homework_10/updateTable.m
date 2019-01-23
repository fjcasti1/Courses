function [Tu,Tv,TY,TR]=updateTable(put,pvt,pYt,Rt,i,Tu,Tv,TY,TR)
r=2;
Fsec=1.25;

% GCI for u
pu=log(abs(put(i-2)-put(i-1))/abs(put(i-1)-put(i)))/log(r);
u_h0=put(i)+(put(i)-put(i-1))/(r^pu-1);
uGCI12=Fsec*abs(1-put(i-1)/put(i))/(r^pu-1)*100;
uGCI23=Fsec*abs(1-put(i-2)/put(i-1))/(r^pu-1)*100;
ubeta=uGCI12*r^pu/uGCI23;
% Include results in a table
Tu(i,3)=pu;
Tu(i,4)=u_h0;
Tu(i,5)=uGCI12;
Tu(i,6)=uGCI23;
Tu(i,7)=ubeta;

% GCI for v
pv=log(abs(pvt(i-2)-pvt(i-1))/abs(pvt(i-1)-pvt(i)))/log(r);
v_h0=pvt(i)+(pvt(i)-pvt(i-1))/(r^pv-1);
vGCI12=Fsec*abs(1-pvt(i-1)/pvt(i))/(r^pv-1)*100;
vGCI23=Fsec*abs(1-pvt(i-2)/pvt(i-1))/(r^pv-1)*100;
vbeta=vGCI12*r^pv/vGCI23;
% Include results in a table
Tv(i,3)=pv;
Tv(i,4)=v_h0;
Tv(i,5)=vGCI12;
Tv(i,6)=vGCI23;
Tv(i,7)=vbeta;

% GCI for Y
pY=log(abs(pYt(i-2)-pYt(i-1))/abs(pYt(i-1)-pYt(i)))/log(r);
Y_h0=pYt(i)+(pYt(i)-pYt(i-1))/(r^pY-1);
YGCI12=Fsec*abs(1-pYt(i-1)/pYt(i))/(r^pY-1)*100;
YGCI23=Fsec*abs(1-pYt(i-2)/pYt(i-1))/(r^pY-1)*100;
Ybeta=YGCI12*r^pY/YGCI23;
% Include results in a table
TY(i,3)=pY;
TY(i,4)=Y_h0;
TY(i,5)=YGCI12;
TY(i,6)=YGCI23;
TY(i,7)=Ybeta;

% GCI for R
pR=log(abs(Rt(i-2)-Rt(i-1))/abs(Rt(i-1)-Rt(i)))/log(r);
R_h0=Rt(i)+(Rt(i)-Rt(i-1))/(r^pR-1);
RGCI12=Fsec*abs(1-Rt(i-1)/Rt(i))/(r^pR-1)*100;
RGCI23=Fsec*abs(1-Rt(i-2)/Rt(i-1))/(r^pR-1)*100;
Rbeta=RGCI12*r^pR/RGCI23;
% Include results in a table
TR(i,3)=pR;
TR(i,4)=R_h0;
TR(i,5)=RGCI12;
TR(i,6)=RGCI23;
TR(i,7)=Rbeta;
end