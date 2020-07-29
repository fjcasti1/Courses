function u = initialCondition(x,M)
u=zeros(M+4,1);
for i=3:M+2
    if (x(i)<4.5 || x(i)>5.5)
        u(i)=0.25+0.5*sin(pi/4*(x(i)-3));
    elseif (x(i)>=4.5 && x(i)<=5.5)
        u(i)=0.25+0.5*sin(pi/4*(x(i)-3))+(1+cos(2*pi*x(i)))*cos(8*pi*x(i));
    end
end
% Periodic Boundary Conditions
u(1)=u(M+1);
u(2)=u(M+2);
u(M+3)=u(3);
u(M+4)=u(4);
end
