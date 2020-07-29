function phi = applyBCs(phi,M,N,opt)
global yu;
global xv;
global xY;
global yY;
switch opt
    case 'u'
        %---------- Boundary Conditions for u----------%
        % Left and right
        for j=1:N+2
            if (yu(j)>=0.5 && yu(j)<=1) % Inlet 1
                phi(1,j)=-24*yu(j)^2+36*yu(j)-12;
            elseif (yu(j)>=1 && yu(j)<=1.5) % Inlet 2
                phi(M+1,j)=24*yu(j)^2-60*yu(j)+36;
            end
        end
        % Top and bottom
        phi(:,1)=-phi(:,2); 
        phi(:,N+2)=-phi(:,N+1);
	case 'v'
        %---------- Boundary Conditions ----------%
        % Top and bottom
        for i=1:M+2
            if (xv(i)>=0.5 && xv(i)<=1) % Inlet 3
                phi(i,N+1)=24*xv(i)^2-36*xv(i)+12;
            elseif (xv(i)>=1.5 && xv(i)<=2.5) % Outlet
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
            if (yY(j)>=0.5 && yY(j)<=1) % Inlet 1
                phi(3,j)=2-phi(4,j);
                phi(2,j)=2-phi(5,j);
                phi(1,j)=2-phi(6,j);
            elseif (yY(j)>=1 && yY(j)<=1.5) % Inlet 2
                phi(M+4,j)=2-phi(M+3,j);
                phi(M+5,j)=2-phi(M+2,j);
                phi(M+6,j)=2-phi(M+1,j);
            end
        end
        %Inlet 3
        for i=3:M+3
            if (xY(i)>=0.5 && xY(i)<=1) % Inlet 3
                phi(i,N+4)=-phi(i,N+3);
                phi(i,N+5)=-phi(i,N+2);
                phi(i,N+6)=-phi(i,N+1);
            end
        end
end
end