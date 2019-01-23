function v=outletCorrection(u,v)
global hx; global hy;
global xv;
global ao; global bo;
s=-sum(u(1,:))*hy+sum(u(end,:))*hy-sum(v(:,1))*hx+sum(v(:,end))*hx;
num=find(xv<=bo,1,'last')-find(xv>=ao,1)+1;
vcorr=s/(num*hx);
for i=find(xv>=ao,1):find(xv<=bo,1,'last')
    v(i,1)=v(i,1)+vcorr;
end
end