function phi = applyBCs(phi,M,N,time,opt)

global yu; 
global xv;
global xY; global yY;
global a1; global b1; global uavg1; global startTime1;
global a2; global b2; global uavg2; global startTime2;
global a3; global b3; global uavg3; global startTime3;
global ao; global bo;

Yin1=1;
Yin2=1;
Yin3=0;
% [Yin1,Yin2,Yin3]=Yin(time);
[c11,c12,c13]=inletCoefficients(uavg1,a1,b1,time,startTime1);% Coefficients for inlet 1
[c21,c22,c23]=inletCoefficients(uavg2,a2,b2,time,startTime2);% Coefficients for inlet 2
[c31,c32,c33]=inletCoefficients(uavg3,a3,b3,time,startTime3);% Coefficients for inlet 3

switch opt
    case 'u'
        %---------- Boundary Conditions for u----------%
        % Left and right
        for j=1:N+2
            if (yu(j)>=a1 && yu(j)<=b1) % Inlet 1
                phi(1,j)=c11*yu(j)^2+c12*yu(j)+c13;
            end
            if (yu(j)>=a2 && yu(j)<=b2) % Inlet 2
                phi(M+1,j)=c21*yu(j)^2+c22*yu(j)+c23;
            end
        end
        % Top and bottom
        phi(:,1)=-phi(:,2); 
        phi(:,N+2)=-phi(:,N+1);
	case 'v'
        %---------- Boundary Conditions ----------%
        % Top and bottom
        for i=1:M+2
            if (xv(i)>=a3 && xv(i)<=b3) % Inlet 3
                phi(i,N+1)=c31*xv(i)^2+c32*xv(i)+c33;
            end
            if (xv(i)>=ao && xv(i)<=bo) % Outlet
                phi(i,1)=(4*phi(i,2)-phi(i,3))/3;
            end
        end
        % Left and right
        phi(1,:)=-phi(2,:); 
        phi(M+2,:)=-phi(M+1,:);
    case 'Y'
        %---------- Boundary Conditions ----------%
        % Zero Neumann for walls and outlet
        %Left
        phi(3,4:N+3)=phi(4,4:N+3); 
        phi(2,4:N+3)=phi(5,4:N+3); 
        phi(1,4:N+3)=phi(6,4:N+3); 
        %Right
        phi(M+4,4:N+3)=phi(M+3,4:N+3); 
        phi(M+5,4:N+3)=phi(M+2,4:N+3); 
        phi(M+6,4:N+3)=phi(M+1,4:N+3); 
        %Bottom
        phi(4:M+3,3)=phi(4:M+3,4);
        phi(4:M+3,2)=phi(4:M+3,5);
        phi(4:M+3,1)=phi(4:M+3,6);
        %Top
        phi(4:M+3,N+4)=phi(4:M+3,N+3);
        phi(4:M+3,N+5)=phi(4:M+3,N+2);
        phi(4:M+3,N+6)=phi(4:M+3,N+1);
        % Dirichlet for the Inlets
        %Inlets 1 and 2
        for j=3:N+3
            if (yY(j)>=a1 && yY(j)<=b1) % Inlet 1
                phi(3,j)=2*Yin1-phi(4,j);
                phi(2,j)=2*Yin1-phi(5,j);
                phi(1,j)=2*Yin1-phi(6,j);
            end
            if (yY(j)>=a2 && yY(j)<=b2) % Inlet 2
                phi(M+4,j)=2*Yin2-phi(M+3,j);
                phi(M+5,j)=2*Yin2-phi(M+2,j);
                phi(M+6,j)=2*Yin2-phi(M+1,j);
            end
        end
        %Inlet 3
        for i=3:M+3
            if (xY(i)>=a3 && xY(i)<=b3) % Inlet 3
                phi(i,N+4)=2*Yin3-phi(i,N+3);
                phi(i,N+5)=2*Yin3-phi(i,N+2);
                phi(i,N+6)=2*Yin3-phi(i,N+1);
            end
        end
end
end