function [u,v,Hu,Hv]=solveBurgers2D(u,v,Hu,Hv,M,N,hx,hy,dt,Re,time)
Huold=Hu;
Hvold=Hv;
[Hu,Hv]=HyperbolicBurgers2D(u,v,M,N,hx,hy);
if time==0 %FCTS for the first time step
    Huold=Hu;
    Hvold=Hv;
end
% % % % % for i=1:M+1
% % % % %     for j=1:N+2
% % % % %         HuoldVector((j-1)*(M+1)+i,1)=Huold(i,j);
% % % % %     end
% % % % % end
% % % % % for i=1:M+2
% % % % %     for j=1:N+1
% % % % %         HvoldVector((j-1)*(M+2)+i,1)=Hvold(i,j);
% % % % %     end
% % % % % end
u=ADI_u(u,M,N,dt,hx,hy,Re,1.5*Hu-0.5*Huold);
v=ADI_v(v,M,N,dt,hx,hy,Re,1.5*Hv-0.5*Hvold);

% % % % % for i=1:M+1
% % % % %         for j=1:N+2
% % % % %             uVector((j-1)*(M+1)+i,1)=u(i,j);
% % % % %         end
% % % % % end
% % % % % for i=1:M+2
% % % % %         for j=1:N+1
% % % % %             vVector((j-1)*(M+2)+i,1)=v(i,j);
% % % % %         end
% % % % % end
% % % % % keyboard
end