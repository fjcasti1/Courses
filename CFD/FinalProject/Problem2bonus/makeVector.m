function [uVector,vVector,YVector,HuVector,HvVector]=makeVector(u,v,Y,Hu,Hv,M,N)
uVector=zeros((M+1)*(N+2),1);
vVector=zeros((M+2)*(N+1),1);
YVector=zeros((M+6)*(N+6),1);
HuVector=zeros((M+1)*(N+2),1);
HvVector=zeros((M+2)*(N+1),1);
for i=1:M+1
    for j=1:N+2
        uVector((j-1)*(M+1)+i,1)=u(i,j);
    end
end
for i=1:M+2
    for j=1:N+1
        vVector((j-1)*(M+2)+i,1)=v(i,j);
    end
end
for i=1:M+6
    for j=1:N+6
        YVector((j-1)*(M+6)+i,1)=Y(i,j);
	end
end
for i=1:M+1
    for j=1:N+2
        HuVector((j-1)*(M+1)+i,1)=Hu(i,j);
    end
end
for i=1:M+2
    for j=1:N+1
        HvVector((j-1)*(M+2)+i,1)=Hv(i,j);
    end
end
end