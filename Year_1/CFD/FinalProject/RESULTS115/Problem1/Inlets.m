clear all
format long
syms varA varB varC
vars=[varA varB varC];
%% Inlet 1
a=0.5;
b=1;
uavg=1;
eqns1=[(b^3-a^3)*varA/3+(b^2-a^2)*varB/2+(b-a)*varC==(b-a)*uavg,...
    a^2*varA+a*varB+varC==0, b^2*varA+b*varB+varC==0];
[A1,B1,C1]=solve(eqns1,vars);
A1=double(A1)
B1=double(B1)
C1=double(C1)
y=linspace(a,b,50e5);
x=A1*y.^2+B1*y+C1;
mean(x)
pause
%% Inlet 2
a=1;
b=1.5;
uavg=-1;
eqns2=[(b^3-a^3)*varA/3+(b^2-a^2)*varB/2+(b-a)*varC==(b-a)*uavg,...
    a^2*varA+a*varB+varC==0, b^2*varA+b*varB+varC==0];
[A2,B2,C2]=solve(eqns2,vars);
A2=double(A2)
B2=double(B2)
C2=double(C2)
x=linspace(a,b,50e5);
y=A2*x.^2+B2*x+C2;
mean(y)
pause
%% Inlet 3
a=0.5;
b=1;
vavg=-1;
eqns3=[(b^3-a^3)*varA/3+(b^2-a^2)*varB/2+(b-a)*varC==(b-a)*vavg,...
    a^2*varA+a*varB+varC==0, b^2*varA+b*varB+varC==0];
[A3,B3,C3]=solve(eqns3,vars);
A3=double(A3)
B3=double(B3)
C3=double(C3)
x=linspace(a,b,50e5);
y=A3*x.^2+B3*x+C3;
mean(y)
pause