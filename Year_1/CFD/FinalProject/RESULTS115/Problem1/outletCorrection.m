function v=outletCorrection(u,v,xv,hx,hy)
s=-sum(u(1,:))*hy+sum(u(end,:))*hy-sum(v(:,1))*hx+sum(v(:,end))*hx;
num=find(xv<=2.5,1,'last')-find(xv>=1.5,1)+1;
vcorr=s/(num*hx);
for i=find(xv>=1.5,1):find(xv<=2.5,1,'last')
    v(i,1)=v(i,1)+vcorr;
end
end