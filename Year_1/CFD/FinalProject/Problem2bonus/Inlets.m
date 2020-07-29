clear all
format long
syms varA varB varC
vars=[varA varB varC];

IOdata=importdata('IOdata.txt',',');
a1=IOdata.data(1,1);
b1=IOdata.data(1,2);
a2=IOdata.data(2,1);
b2=IOdata.data(2,2);
a3=IOdata.data(3,1);
b3=IOdata.data(3,2);
a4=IOdata.data(4,1);
b4=IOdata.data(4,2);
ao=IOdata.data(5,1);
bo=IOdata.data(5,2);

%% Inlet 1
uavg=1;
eqns1=[(b1^3-a1^3)*varA/3+(b1^2-a1^2)*varB/2+(b1-a1)*varC==(b1-a1)*uavg,...
    a1^2*varA+a1*varB+varC==0, b1^2*varA+b1*varB+varC==0];
[A1,B1,C1]=solve(eqns1,vars);
A1=double(A1)
B1=double(B1)
C1=double(C1)
y=linspace(a1,b1,50e5);
x=A1*y.^2+B1*y+C1;
mean(x)
pause
%% Inlet 2
uavg=-1;
eqns2=[(b2^3-a2^3)*varA/3+(b2^2-a2^2)*varB/2+(b2-a2)*varC==(b2-a2)*uavg,...
    a2^2*varA+a2*varB+varC==0, b2^2*varA+b2*varB+varC==0];
[A2,B2,C2]=solve(eqns2,vars);
A2=double(A2)
B2=double(B2)
C2=double(C2)
x=linspace(a2,b2,50e5);
y=A2*x.^2+B2*x+C2;
mean(y)
pause
%% Inlet 3
vavg=-1;
eqns3=[(b3^3-a3^3)*varA/3+(b3^2-a3^2)*varB/2+(b3-a3)*varC==(b3-a3)*vavg,...
    a3^2*varA+a3*varB+varC==0, b3^2*varA+b3*varB+varC==0];
[A3,B3,C3]=solve(eqns3,vars);
A3=double(A3)
B3=double(B3)
C3=double(C3)
x=linspace(a3,b3,50e5);
y=A3*x.^2+B3*x+C3;
mean(y)
pause
%% Inlet 4
vavg=-1;
eqns3=[(b4^3-a4^3)*varA/3+(b4^2-a4^2)*varB/2+(b4-a4)*varC==(b4-a4)*vavg,...
    a4^2*varA+a4*varB+varC==0, b4^2*varA+b4*varB+varC==0];
[A4,B4,C4]=solve(eqns3,vars);
A4=double(A4)
B4=double(B4)
C4=double(C4)
x=linspace(a4,b4,50e5);
y=A4*x.^2+B4*x+C4;
mean(y)
pause