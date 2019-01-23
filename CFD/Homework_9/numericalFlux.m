function h = numericalFlux(E,Phi,M)
h=zeros(M+1,1);
for i=2:M+2
    h(i-1)=0.5*(E(i+1)+E(i)+Phi(i-1));
end