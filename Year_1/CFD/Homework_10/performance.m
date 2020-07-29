function R=performance(Y,M,N,hx,hy,L,H)
R=sum(sum(Y(4:M+3,4:N+3).*(1-Y(4:M+3,4:N+3))))*hx*hy/(H*L);
end